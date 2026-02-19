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
    def __init__(self, container_name='octave_runner', verbose=False):
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

        # Shared directory for temporary files
        self.shared_dir = os.path.join(self.project_root, 'docker_setup', 'octave_shared')
        os.makedirs(self.shared_dir, exist_ok=True)

        self.script_path = os.path.join(self.shared_dir, 'run.m')
        self.result_path = os.path.join(self.shared_dir, 'result.mat')
        self.commands = []  # Temporary commands cleared after each feval
        self.persistent_commands = []  # Persistent commands (like addpath) that stay
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

        # Check if the container exists
        try:
            self.container = self.client.containers.get(container_name)
            if self.container.status != 'running':
                print(f"Starting existing Docker container '{container_name}'...")
                self.container.start()
            else:
                print(f"Using existing Docker container '{container_name}'")
        except docker.errors.NotFound:
            print(f"Creating new Docker container '{container_name}' with Octave...")
            # Mount the entire project directory to /workspace in the container
            self.container = self.client.containers.run(
                'basisremy-octave:latest',
                name=container_name,
                command='tail -f /dev/null',
                volumes={self.project_root: {'bind': '/workspace', 'mode': 'rw'}},
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

        try:
            # Try to get existing image
            self.client.images.get(image_name)
            print(f"✓ Using existing Docker image '{image_name}'")
        except docker.errors.ImageNotFound:
            # Image doesn't exist, build it
            print(f"Building BasisREMY Octave Docker image (this may take a few minutes)...")
            print("=" * 80)
            try:
                # Get the docker_setup directory path
                dockerfile_dir = os.path.dirname(os.path.abspath(__file__))

                # Build the image from the Dockerfile with streaming output
                build_logs = self.client.api.build(
                    path=os.path.dirname(dockerfile_dir),  # Project root
                    dockerfile=os.path.join('docker_setup', 'dockerfile'),
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
                print(f"✓ BasisREMY Octave Docker image built successfully")
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
        # Normalize path - remove leading './'
        normalized_path = path.replace('\\', '/').lstrip('./')
        # Return genpath expression - this will be evaluated in Octave
        return f"genpath('{normalized_path}')"

    def addpath(self, path_or_genpath_result):
        """Add a path to Octave's search path (persistent)."""
        if isinstance(path_or_genpath_result, str):
            if 'genpath(' in path_or_genpath_result:
                # This is a genpath result, use it directly
                self.persistent_commands.append(f"addpath({path_or_genpath_result});")
            else:
                # This is a regular path - normalize it
                normalized_path = path_or_genpath_result.replace('\\', '/').lstrip('./')
                self.persistent_commands.append(f"addpath('{normalized_path}');")

    def set_verbose(self, verbose):
        """Enable or disable verbose output."""
        self.verbose = bool(verbose)
        if self.verbose:
            print("✓ Docker Octave verbose mode enabled")

    def check_running_processes(self):
        """Check for existing Octave processes in the container."""
        try:
            result = self.container.exec_run("pgrep -a octave-cli")
            if result.exit_code == 0:
                processes = result.output.decode().strip().split('\n')
                return [p for p in processes if p]
            return []
        except Exception:
            return []

    def kill_running_processes(self):
        """Kill all running Octave processes in the container."""
        try:
            result = self.container.exec_run("pkill -9 octave-cli")
            if self.verbose:
                print("✓ Killed existing Octave processes")
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
                assigns.append(f"{var} = '{normalized_arg}';")
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
                    cell_items = ', '.join(f"'{item}'" for item in arg)
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

        # Save to shared directory (relative to /workspace)
        result_file_rel = os.path.relpath(self.result_path, self.project_root).replace('\\', '/')
        save = f"save('-v7', '{result_file_rel}', {store_vars_str});"

        # Build the complete script - include persistent commands first
        code = '\n'.join(self.persistent_commands + assigns + self.commands + [call, save])

        if show_output:
            print(f"\nGenerated Octave script:")
            print(f"{'-'*80}")
            # Show only the key parts if verbose
            print('\n'.join(self.persistent_commands[:3]) + '\n...')
            print(call)
            print(save)
            print(f"{'-'*80}")

        # Write script to shared directory
        with open(self.script_path, 'w') as f:
            f.write(code)

        if show_output:
            print(f"\n✓ Script written to: {self.script_path}")
            print(f"⏳ Executing Octave in Docker container...")

        # Execute in container - script path relative to /workspace
        script_rel = os.path.relpath(self.script_path, self.project_root).replace('\\', '/')

        if show_output:
            print(f"   Command: octave-cli {script_rel}")
            print(f"{'-'*80}")

        # Check for existing Octave processes
        existing_check = self.container.exec_run("pgrep octave-cli")
        if existing_check.exit_code == 0:
            existing_pids = existing_check.output.decode().strip().split('\n')
            if existing_pids and existing_pids[0]:
                print(f"⚠️  Warning: Found {len(existing_pids)} existing Octave process(es) running!")
                print(f"   PIDs: {', '.join(existing_pids)}")
                print(f"   This may slow down your simulation significantly.")
                print(f"   Consider killing them with: docker exec octave_runner pkill -9 octave-cli")
                print(f"{'-'*80}")

        # Add helpful message for long-running simulations
        if show_output:
            print("📝 Note: Basis set simulations can take several minutes per metabolite.")
            print("   The process is running if you see this message - please be patient!")
            print(f"{'-'*80}")

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
            print(f"✓ Octave execution completed successfully")
            print(f"⏳ Loading results from {result_file_rel}...")

        # Load results
        try:
            # Use squeeze_me=True to remove singleton dimensions from arrays
            # Use struct_as_record=False to get more intuitive struct access
            mat = scipy.io.loadmat(self.result_path, squeeze_me=True, struct_as_record=False)
            if show_output:
                print(f"✓ Results loaded successfully")
                print(f"{'='*80}\n")
        except Exception as e:
            print(f"✗ Failed to load results: {e}")
            raise RuntimeError(f"Failed to load Octave results: {e}")

        # Clear commands for next execution
        self.commands = []

        # Return results
        if store_as:
            return mat[store_as]

        if nout == 1:
            return mat[result_vars[0]]
        else:
            return tuple(mat[v] for v in result_vars)

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


