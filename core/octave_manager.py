####################################################################################################
#                                         octave_manager.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 11/02/26                                                                                #
#                                                                                                  #
# Purpose: Manages Octave runtime environment with Docker and local fallback support.              #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import os
import subprocess
import shutil


#**************************************************************************************************#
#                                         OctaveManager                                            #
#**************************************************************************************************#
#                                                                                                  #
# Manages the Octave runtime environment. Attempts to use Docker first, then falls back to local   #
# Octave installation if Docker is not available.                                                  #
#                                                                                                  #
#**************************************************************************************************#
class OctaveManager:
    """
    Manages Octave runtime with Docker-first approach and local fallback.

    Priority order:
    1. Try Docker (if available and running)
    2. Fall back to local Octave installation
    3. Raise error with helpful instructions if neither is available
    """

    def __init__(self, verbose=False):
        self.runtime_type = None  # Will be 'docker', 'local', or None
        self.octave_instance = None
        self.docker_available = False
        self.local_octave_available = False

        # Check for environment variable to enable verbose mode
        env_verbose = os.environ.get('BASISREMY_VERBOSE', '').lower() in ('1', 'true', 'yes')
        self.verbose = bool(verbose) or env_verbose

        if self.verbose and env_verbose:
            print("✓ Verbose mode enabled via BASISREMY_VERBOSE environment variable")

    def check_docker_availability(self):
        """Check if Docker is installed and running."""
        try:
            import docker

            # Try to connect to Docker - handle different socket locations
            # First try default (works for most setups)
            try:
                client = docker.from_env()
                client.ping()
                self.docker_available = True
                return True
            except Exception:
                # Try OrbStack socket location (macOS)
                try:
                    orbstack_socket = os.path.expanduser('~/.orbstack/run/docker.sock')
                    if os.path.exists(orbstack_socket):
                        client = docker.DockerClient(base_url=f'unix://{orbstack_socket}')
                        client.ping()
                        self.docker_available = True
                        return True
                except Exception:
                    pass

                # Try other common locations
                for socket_path in ['/var/run/docker.sock',
                                   os.path.expanduser('~/Library/Containers/com.docker.docker/Data/docker.sock')]:
                    try:
                        if os.path.exists(socket_path):
                            client = docker.DockerClient(base_url=f'unix://{socket_path}')
                            client.ping()
                            self.docker_available = True
                            return True
                    except Exception:
                        continue

                return False

        except ImportError:
            # Docker Python package not installed
            return False
        except Exception as e:
            # Docker not running or other error
            return False

    def check_local_octave_availability(self):
        """Check if Octave is installed locally."""
        try:
            # Check if octave or octave-cli is in PATH
            octave_cmd = shutil.which('octave-cli') or shutil.which('octave')
            if octave_cmd:
                # Try to run a simple command to verify it works
                result = subprocess.run(
                    [octave_cmd, '--eval', 'disp("test")'],
                    capture_output=True,
                    timeout=5,
                    text=True
                )
                if result.returncode == 0:
                    self.local_octave_available = True
                    return True
        except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
            pass
        return False

    def initialize_octave(self, prefer_docker=True, verbose=None):
        """
        Initialize Octave runtime with fallback logic.

        Args:
            prefer_docker: If True, try Docker first. If False, try local first.
            verbose: If provided, override the manager's verbose flag for this call.

        Returns:
            octave_instance: Either DockerOctave or Oct2Py instance

        Raises:
            RuntimeError: If neither Docker nor local Octave is available
        """
        # Allow per-call verbose override
        if verbose is not None:
            self.verbose = bool(verbose)

        if prefer_docker:
            # Try Docker first
            if self.check_docker_availability():
                return self._initialize_docker()
            # Fall back to local
            elif self.check_local_octave_availability():
                return self._initialize_local()
            else:
                raise RuntimeError(self._get_installation_instructions())
        else:
            # Try local first
            if self.check_local_octave_availability():
                return self._initialize_local()
            # Fall back to Docker
            elif self.check_docker_availability():
                return self._initialize_docker()
            else:
                raise RuntimeError(self._get_installation_instructions())

    def _initialize_docker(self):
        """Initialize Docker-based Octave."""
        try:
            from docker_setup.docker_octave import DockerOctave

            if self.verbose:
                print("Initializing Docker-based Octave runtime (verbose mode)...")

            # Pass verbose flag to DockerOctave
            self.octave_instance = DockerOctave(verbose=self.verbose)
            self.runtime_type = 'docker'

            if self.verbose:
                print("✓ Docker Octave runtime initialized successfully")
            else:
                print("✓ Using Docker for Octave runtime")

            return self.octave_instance
        except Exception as e:
            raise RuntimeError(f"Failed to initialize Docker Octave: {e}")

    def _initialize_local(self):
        """Initialize local Octave installation."""
        try:
            from oct2py import Oct2Py
            self.octave_instance = Oct2Py()
            self.runtime_type = 'local'
            print("✓ Using local Octave installation")
            return self.octave_instance
        except Exception as e:
            raise RuntimeError(f"Failed to initialize local Octave: {e}")

    def _get_installation_instructions(self):
        """Generate helpful installation instructions based on what's missing."""
        instructions = [
            "\n" + "="*80,
            "Octave Runtime Not Available",
            "="*80,
            "",
            "BasisREMY requires Octave to run simulations. You have two options:",
            "",
        ]

        # Check what's available and provide specific instructions
        docker_checked = self.check_docker_availability()
        octave_checked = self.check_local_octave_availability()

        instructions.append("OPTION 1: Use Docker (Recommended)")
        instructions.append("-" * 40)
        if not docker_checked:
            instructions.extend([
                "Docker is not installed or not running.",
                "",
                "To install Docker:",
                "  • macOS: Download from https://www.docker.com/products/docker-desktop",
                "  • Linux: sudo apt-get install docker.io (Ubuntu/Debian)",
                "           sudo yum install docker (RedHat/CentOS)",
                "  • Windows: Download from https://www.docker.com/products/docker-desktop",
                "",
                "After installation, make sure Docker is running and try again.",
                "You also need to install the Python docker package:",
                "  pip install docker",
                "",
            ])
        else:
            instructions.extend([
                "✓ Docker is available!",
                "",
            ])

        instructions.append("OPTION 2: Install Octave Locally")
        instructions.append("-" * 40)
        if not octave_checked:
            instructions.extend([
                "Local Octave installation not found.",
                "",
                "To install Octave:",
                "  • macOS: brew install octave",
                "  • Linux: sudo apt-get install octave (Ubuntu/Debian)",
                "           sudo yum install octave (RedHat/CentOS)",
                "  • Windows: Download from https://www.gnu.org/software/octave/",
                "",
                "You also need to install the Python oct2py package:",
                "  pip install oct2py",
                "",
            ])
        else:
            instructions.extend([
                "✓ Local Octave is available!",
                "",
            ])

        instructions.append("="*80)

        return "\n".join(instructions)

    def get_runtime_info(self):
        """Get information about the current runtime."""
        return {
            'runtime_type': self.runtime_type,
            'docker_available': self.docker_available,
            'local_octave_available': self.local_octave_available,
            'verbose': self.verbose,
        }
