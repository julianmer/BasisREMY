####################################################################################################
#                                          docker_octave.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: A. Wright (andrew.wright@utsouthwestern.edu)                                            #
#          J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 14/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines a Docker-based interface to run Octave commands.                                #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import docker
import glob
import itertools
import numpy as np
import os
import scipy.io


#**************************************************************************************************#
#                                           DockerOctave                                           #
#**************************************************************************************************#
#                                                                                                  #
# A class to run Octave commands inside a Docker container. It manages the container lifecycle,    #
# script generation, and result retrieval.                                                         #
#                                                                                                  #
#**************************************************************************************************#
class DockerOctave:
    # must match the LABEL in docker/dockerfile; bump both to force a rebuild
    IMAGE_VERSION = 'lean-2'

    # Default ``addpath`` prefix for the bundled adapter scripts, relative to
    # ``/workspace`` (the container working dir). ``__init__`` overrides this
    # per instance based on where the adapters sit relative to the mounted
    # project root (see the volume-mount block below); backends read the
    # effective value via :func:`basisremy.core.paths.octave_adapters_base`.
    ADAPTERS_MOUNT = 'adapters'

    def __init__(self, container_name=None, verbose=False):
        """
        Initialize Docker-based Octave runtime.

        Args:
            container_name: Name for the Docker container
            verbose: Enable verbose logging for debugging

        Mounts the entire current working directory to /workspace in the container
        so that all file paths work transparently.
        """
        self.verbose = bool(verbose)

        # Get the project root directory (where we're running from)
        self.project_root = os.getcwd()

        # Scratch directory for the generated run.m / result.mat. It MUST live
        # under the working directory because that is what gets mounted into the
        # container at /workspace (see the volume mount below).
        self.shared_dir = os.path.join(self.project_root, '.octave_shared')
        os.makedirs(self.shared_dir, exist_ok=True)

        # One script/result pair per process AND per call: two BasisREMY
        # processes in the same project folder must not overwrite each other's
        # run, and neither must the threads of one process simulating several
        # metabolites at once (see Backend.simulate_in_parallel). The names
        # share the prefix `run_<pid>_`, which is what the own-process checks
        # and kills match on.
        self._scratch_prefix = f'{os.getpid()}_'
        self._call_counter = itertools.count()
        self.script_path = os.path.join(self.shared_dir, f'run_{self._scratch_prefix}0.m')
        self.result_path = os.path.join(self.shared_dir, f'result_{self._scratch_prefix}0.mat')
        import atexit
        atexit.register(self._remove_scratch_files)
        self.commands = []  # Temporary commands cleared after each feval
        self.persistent_commands = []  # Persistent commands (like addpath) that stay
        if container_name is None:
            # Per-project container: sessions running from different working
            # directories mount different trees and would otherwise destroy
            # each other's container when refreshing stale volume mounts.
            import hashlib
            # sha256, not md5: md5 raises on FIPS-enabled hosts
            digest = hashlib.sha256(self.project_root.encode()).hexdigest()[:8]
            container_name = f'octave_runner_{digest}'
        self.container_name = container_name

        # Try to connect to Docker - handle different socket locations
        try:
            # First try default connection
            self.client = docker.from_env()
            self.client.ping()  # Verify connection works
        except Exception:
            # Try OrbStack socket location (macOS)
            try:
                orbstack_socket = os.path.expanduser('~/.orbstack/run/docker.sock')
                if os.path.exists(orbstack_socket):
                    self.client = docker.DockerClient(base_url=f'unix://{orbstack_socket}')
                    self.client.ping()
                else:
                    raise Exception("OrbStack socket not found")
            except Exception:
                # Try other common socket locations
                socket_found = False
                for socket_path in ['/var/run/docker.sock',
                                   os.path.expanduser('~/Library/Containers/com.docker.docker/Data/docker.sock')]:
                    try:
                        if os.path.exists(socket_path):
                            self.client = docker.DockerClient(base_url=f'unix://{socket_path}')
                            self.client.ping()
                            socket_found = True
                            break
                    except Exception:
                        continue

                if not socket_found:
                    raise RuntimeError(
                        "Failed to connect to Docker. Please ensure Docker is installed and running.\n"
                        "Tried locations: docker.from_env(), OrbStack, /var/run/docker.sock"
                    )

        # Pull Octave image if not present
        self._ensure_octave_image()

        # The bundled adapter scripts live inside the installed package. When
        # the package sits inside the mounted project directory (the usual
        # source-checkout case) they are already visible through the single
        # /workspace mount, so we must NOT add a second bind-mount underneath
        # /workspace: Docker Desktop for Mac silently shadows such nested
        # mounts, leaving /workspace/adapters empty and the adapter scripts
        # unreachable. Only when the adapters live outside the project root do
        # we mount them separately.
        from basisremy.core.paths import ADAPTERS_DIR

        volumes = {
            self.project_root: {'bind': '/workspace', 'mode': 'rw'},
        }
        adapters_rel = os.path.relpath(str(ADAPTERS_DIR), self.project_root)
        if not adapters_rel.startswith(os.pardir) and not os.path.isabs(adapters_rel):
            # Inside the project: reach the adapters via the /workspace mount.
            self.ADAPTERS_MOUNT = adapters_rel.replace(os.sep, '/')
        else:
            # Outside the project: bind-mount them in at a dedicated location.
            volumes[str(ADAPTERS_DIR)] = {'bind': '/workspace/adapters', 'mode': 'ro'}
            self.ADAPTERS_MOUNT = 'adapters'

        # Check if the container exists
        try:
            self.container = self.client.containers.get(container_name)
            # A container created by an older version — or one whose bind-mounts
            # point at a stale host path — would not expose the current scripts.
            # Validate every expected mount against the running container and
            # recreate on any mismatch. Crucially, this also tears down legacy
            # containers that still carry the broken nested /workspace/adapters
            # mount (silently shadowed on Docker Desktop for Mac).
            mounts = self.container.attrs.get('Mounts', [])
            actual = {
                m.get('Destination'): os.path.realpath(m.get('Source', ''))
                for m in mounts
            }
            expected = {
                spec['bind']: os.path.realpath(src)
                for src, spec in volumes.items()
            }
            mounts_ok = actual == expected
            # a container created from a previous image keeps running the old
            # Octave — recreate it on the current image as well
            current_image = self.client.images.get('basisremy-octave:latest').id
            image_ok = self.container.image.id == current_image
            if not (mounts_ok and image_ok):
                print(
                    f"Recreating Docker container '{container_name}' to refresh "
                    + ("the Octave image..." if mounts_ok
                       else "stale, missing, or obsolete volume mounts...")
                )
                self.container.remove(force=True)
                raise docker.errors.NotFound('recreate')
            if self.container.status != 'running':
                print(f"Starting existing Docker container '{container_name}'...")
                self.container.start()
            else:
                print(f"Using existing Docker container '{container_name}'")
        except docker.errors.NotFound:
            # Names are per-project now; a fixed-name container from an older
            # version keeps running (holding its mounts) until removed by hand.
            # Never auto-remove it — an old BasisREMY session may still use it.
            try:
                self.client.containers.get('octave_runner')
                print("ℹ️  A legacy 'octave_runner' container from an older "
                      "BasisREMY is still present. Once no old version is "
                      "running, remove it with: docker rm -f octave_runner")
            except Exception:
                pass
            print(f"Creating new Docker container '{container_name}' with Octave...")
            # Mount the project directory (plus the bundled adapters) into the
            # container; working dir is /workspace so relative paths resolve.
            self.container = self.client.containers.run(
                'basisremy-octave:latest',
                name=container_name,
                command='tail -f /dev/null',
                volumes=volumes,
                working_dir='/workspace',
                detach=True
            )
            print(f"✓ Docker container '{container_name}' created successfully")

    def _ensure_octave_image(self):
        """Build or get BasisREMY Octave Docker image"""
        image_name = 'basisremy-octave:latest'

        # TODO: In the future, attempt to pull a prebuilt image from a registry before building locally
        # This will avoid long local builds and provide a better user experience:
        #   try:
        #       print(f"Pulling prebuilt Docker image '{image_name}'...")
        #       self.client.images.pull('yourregistry/basisremy-octave:latest')
        #       print(f"✓ Using prebuilt Docker image '{image_name}'")
        #       return
        #   except docker.errors.ImageNotFound:
        #       print("Prebuilt image not found, building locally...")

        stale_image = None
        try:
            # Try to get existing image
            image = self.client.images.get(image_name)
            if image.labels.get('org.basisremy.octave-image') != self.IMAGE_VERSION:
                # an image built from an older dockerfile (the original one
                # carried the whole Python stack and weighed 16 GB) — rebuild
                print(f"Docker image '{image_name}' comes from an older BasisREMY — "
                      "rebuilding the lean Octave image...")
                stale_image = image
                raise docker.errors.ImageNotFound(image_name)
            print(f"✓ Using existing Docker image '{image_name}'")
        except docker.errors.ImageNotFound:
            # Image doesn't exist, build it
            print("Building BasisREMY Octave Docker image (this may take a few minutes)...")
            print("=" * 80)
            try:
                # The dockerfile ships next to this module inside the package.
                dockerfile_dir = os.path.dirname(os.path.abspath(__file__))

                # Build the image from the Dockerfile with streaming output
                build_logs = self.client.api.build(
                    path=dockerfile_dir,  # Build context: basisremy/docker/
                    dockerfile='dockerfile',
                    tag=image_name,
                    rm=True,  # Remove intermediate containers
                    decode=True  # Decode JSON stream
                )

                # Stream build progress in real-time
                for log in build_logs:
                    if 'stream' in log:
                        print(log['stream'], end='')
                    elif 'error' in log:
                        raise docker.errors.BuildError(log['error'], build_logs)
                    elif 'status' in log:
                        print(log['status'])

                print("=" * 80)
                print("✓ BasisREMY Octave Docker image built successfully")
                if stale_image is not None:
                    # the old image is untagged now; free its space (containers
                    # still using it are recreated by __init__)
                    try:
                        for c in self.client.containers.list(all=True,
                                                             filters={'ancestor': stale_image.id}):
                            c.remove(force=True)
                        self.client.images.remove(stale_image.id, force=True)
                        print("✓ Removed the previous Octave image")
                    except Exception as exc:  # noqa: BLE001
                        print(f"  (previous image left in place: {exc})")
            except docker.errors.BuildError as e:
                print("=" * 80)
                raise RuntimeError(
                    f"Failed to build BasisREMY Octave Docker image.\n"
                    f"Error: {e}"
                )
            except Exception as e:
                print("=" * 80)
                raise RuntimeError(
                    f"Failed to build BasisREMY Octave Docker image.\n"
                    f"Error: {e}"
                )

    def eval(self, cmd):
        """Execute an Octave command (persistent - stays for all feval calls)."""
        self.persistent_commands.append(cmd)

    def genpath(self, path):
        """
        Generate path string including subdirectories.
        Returns the path string that can be used with addpath.
        This mimics Octave's genpath function.
        """
        # Normalize path - remove a leading './' (lstrip would also strip a
        # leading '/' from absolute paths and mangle '../')
        normalized_path = path.replace('\\', '/')
        if normalized_path.startswith('./'):
            normalized_path = normalized_path[2:]
        # Return genpath expression - this will be evaluated in Octave
        return f"genpath('{self._q(normalized_path)}')"

    @staticmethod
    def _q(s: str) -> str:
        """Escape a string for an Octave single-quoted literal ('' = ')."""
        return str(s).replace("'", "''")

    def addpath(self, path_or_genpath_result):
        """Add a path to Octave's search path (persistent)."""
        if isinstance(path_or_genpath_result, str):
            if 'genpath(' in path_or_genpath_result:
                # This is a genpath result, use it directly
                self.persistent_commands.append(f"addpath({path_or_genpath_result});")
            else:
                # This is a regular path - normalize it
                normalized_path = path_or_genpath_result.replace('\\', '/')
                if normalized_path.startswith('./'):
                    normalized_path = normalized_path[2:]
                self.persistent_commands.append(
                    f"addpath('{self._q(normalized_path)}');")

    def set_verbose(self, verbose):
        """Enable or disable verbose output."""
        self.verbose = bool(verbose)
        if self.verbose:
            print("✓ Docker Octave verbose mode enabled")

    def _own_script(self) -> str:
        # prefix shared by every script this process writes (run_<pid>_<n>.m)
        return f'run_{self._scratch_prefix}'

    def check_running_processes(self, own_only=True):
        """Octave processes in the container — by default only the one(s)
        running this instance's script. Other BasisREMY sessions share the
        container and must be left alone."""
        try:
            result = self.container.exec_run("pgrep -a octave-cli")
            if result.exit_code != 0:
                return []
            processes = [p for p in result.output.decode().strip().split('\n') if p]
            if own_only:
                processes = [p for p in processes if self._own_script() in p]
            return processes
        except Exception:
            return []

    def kill_running_processes(self):
        """Kill the Octave process running this instance's script — never the
        simulations of other sessions in the shared container."""
        try:
            self.container.exec_run(f"pkill -9 -f {self._own_script()}")
            if self.verbose:
                print("✓ Killed this session's Octave process")
            return True
        except Exception as e:
            if self.verbose:
                print(f"⚠️  Failed to kill processes: {e}")
            return False

    def feval(self, func_path, *func_args, nout=1, store_as=None, verbose=False, **kwargs):
        """
        Evaluate an Octave function with arguments.

        Args:
            func_path: Name of the Octave function to call
            *func_args: Arguments to pass to the function
            nout: Number of output arguments
            store_as: Variable name to store result as (optional)
            verbose: Print Octave output (overrides instance verbose setting)

        Returns:
            Result(s) from the Octave function
        """
        # Allow per-call verbose override or use instance setting
        show_output = verbose or self.verbose

        if show_output:
            print(f"\n{'='*80}")
            print(f"Docker Octave: Executing {func_path}()")
            print(f"{'='*80}")

        arg_vars = []
        assigns = []

        for i, arg in enumerate(func_args):
            var = f'arg{i}'
            arg_vars.append(var)

            if isinstance(arg, str):
                # Normalize file paths in string arguments
                normalized_arg = arg.replace('\\', '/')
                # Remove leading './' to work from /workspace
                if normalized_arg.startswith('./'):
                    normalized_arg = normalized_arg[2:]
                assigns.append(f"{var} = '{self._q(normalized_arg)}';")
                if show_output:
                    print(f"  arg{i} (str): {normalized_arg}")
            elif isinstance(arg, bool):
                assigns.append(f"{var} = {int(arg)};")
                if show_output:
                    print(f"  arg{i} (bool): {arg}")
            elif isinstance(arg, (int, float)):
                assigns.append(f"{var} = {arg};")
                if show_output:
                    print(f"  arg{i} (num): {arg}")
            elif isinstance(arg, list):
                # Check if it's a list of strings (e.g., metabolite names)
                if arg and isinstance(arg[0], str):
                    # Create Octave cell array
                    cell_items = ', '.join(f"'{self._q(item)}'" for item in arg)
                    assigns.append(f"{var} = {{{cell_items}}};")
                    if show_output:
                        print(f"  arg{i} (list): {len(arg)} items")
                else:
                    # Numeric list - create numeric array
                    assigns.append(f"{var} = [{', '.join(map(str, arg))}];")
                    if show_output:
                        print(f"  arg{i} (array): {len(arg)} elements")
            elif isinstance(arg, np.ndarray):
                assigns.append(f"{var} = [{', '.join(map(str, np.ravel(arg)))}];")
                if show_output:
                    print(f"  arg{i} (ndarray): shape {arg.shape}")
            elif arg is None:
                # Handle None as empty matrix []
                assigns.append(f"{var} = [];")
                if show_output:
                    print(f"  arg{i} (None): []")
            else:
                raise TypeError(f'Unsupported argument type: {type(arg)}')

        # Prepare output variables
        result_vars = [f'result{i}' for i in range(nout)] if isinstance(nout, int) and nout > 1 else ['result']
        call = f"[{', '.join(result_vars)}] = {func_path}({', '.join(arg_vars)});"

        # Determine variables to save
        store_vars = [store_as] if store_as else result_vars
        store_vars_str = ', '.join(repr(v) for v in store_vars)

        # Per-call script / result files (concurrent fevals from several threads
        # must not share them); relative to /workspace inside the container.
        call_id = next(self._call_counter)
        script_path = os.path.join(self.shared_dir, f'run_{self._scratch_prefix}{call_id}.m')
        result_path = os.path.join(self.shared_dir, f'result_{self._scratch_prefix}{call_id}.mat')
        result_file_rel = os.path.relpath(result_path, self.project_root).replace('\\', '/')
        save = f"save('-v7', '{result_file_rel}', {store_vars_str});"

        # Build the complete script - include persistent commands first
        code = '\n'.join(self.persistent_commands + assigns + self.commands + [call, save])

        if show_output:
            print("\nGenerated Octave script:")
            print(f"{'-'*80}")
            # Show only the key parts if verbose
            print('\n'.join(self.persistent_commands[:3]) + '\n...')
            print(call)
            print(save)
            print(f"{'-'*80}")

        # Write script to shared directory
        with open(script_path, 'w') as f:
            f.write(code)

        if show_output:
            print(f"\n✓ Script written to: {script_path}")
            print("⏳ Executing Octave in Docker container...")

        # Execute in container - script path relative to /workspace
        script_rel = os.path.relpath(script_path, self.project_root).replace('\\', '/')

        if show_output:
            print(f"   Command: octave-cli {script_rel}")
            print(f"{'-'*80}")

        # Other BasisREMY sessions may be simulating in this container: say
        # so (both runs get slower) but never touch their processes.
        others = [p for p in self.check_running_processes(own_only=False)
                  if self._own_script() not in p]
        if others:
            print(f"ℹ️  {len(others)} other Octave process(es) running in this "
                  f"container (another BasisREMY session?) — both will be slower.")
            print(f"{'-'*80}")

        # Add helpful message for long-running simulations
        if show_output:
            print("📝 Note: Basis set simulations can take several minutes per metabolite.")
            print("   The process is running if you see this message - please be patient!")
            print(f"{'-'*80}")

        # A result left by the previous run must never pass for this run's
        # output if Octave exits without reaching its save().
        try:
            os.remove(result_path)
        except OSError:
            pass
        exit_code, output = self.container.exec_run(f"octave-cli {script_rel}")

        if show_output or exit_code != 0:
            output_text = output.decode()
            if output_text.strip():
                print("Octave output:")
                print(output_text)
            else:
                if show_output:
                    print("(No output from Octave)")

        if exit_code != 0:
            print(f"{'-'*80}")
            print(f"✗ Octave execution failed with exit code {exit_code}")
            print(f"{'-'*80}")
            raise RuntimeError(f"Octave execution failed with exit code {exit_code}")

        if show_output:
            print(f"{'-'*80}")
            print("✓ Octave execution completed successfully")
            print(f"⏳ Loading results from {result_file_rel}...")

        # Load results
        try:
            # Use squeeze_me=True to remove singleton dimensions from arrays
            # Use struct_as_record=False to get more intuitive struct access
            mat = scipy.io.loadmat(result_path, squeeze_me=True, struct_as_record=False)
            if show_output:
                print("✓ Results loaded successfully")
                print(f"{'='*80}\n")
        except Exception as e:
            print(f"✗ Failed to load results: {e}")
            raise RuntimeError(f"Failed to load Octave results: {e}")
        finally:
            for f in (script_path, result_path):
                try:
                    os.remove(f)
                except OSError:
                    pass

        # Clear commands for next execution
        self.commands = []

        # Return results
        if store_as:
            return mat[store_as]

        if nout == 1:
            return mat[result_vars[0]]
        else:
            return tuple(mat[v] for v in result_vars)

    def _remove_scratch_files(self):
        for pattern in (f'run_{self._scratch_prefix}*.m', f'result_{self._scratch_prefix}*.mat'):
            for p in glob.glob(os.path.join(self.shared_dir, pattern)):
                try:
                    os.remove(p)
                except OSError:
                    pass

    def exit(self):
        """Clear command buffers."""
        self.commands = []
        # Don't clear persistent_commands - they should stay

    def __del__(self):
        """Cleanup when object is destroyed."""
        try:
            if hasattr(self, 'container'):
                # Don't stop the container - it can be reused
                # Just clean up commands
                self.commands = []
        except:
            pass

    def stop_container(self):
        """Stop and remove the Docker container (call this when completely done)."""
        try:
            if hasattr(self, 'container'):
                print(f"Stopping Docker container '{self.container_name}'...")
                self.container.stop()
                self.container.remove()
                print(f"✓ Docker container '{self.container_name}' stopped and removed")
        except Exception as e:
            print(f"Warning: Failed to stop container: {e}")

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - but don't stop container to allow reuse."""
        self.exit()
        return False


