# Octave Runtime Setup Guide

BasisREMY uses Octave to run MRS basis set simulations. The system automatically tries to use the best available option, but you need at least one of the following installed.

## Quick Start

BasisREMY will automatically:
1. **Try Docker first** (recommended for most users)
2. **Fall back to local Octave** if Docker is not available
3. **Show helpful instructions** if neither is available

## Option 1: Docker (Recommended)

Docker is the easiest option and works the same across all platforms.

### Why Docker?
- ✅ No manual Octave installation needed
- ✅ Consistent environment across all platforms
- ✅ Automatic dependency management
- ✅ Easy to set up

### Installing Docker

#### macOS
1. Download [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop)
2. Install and start Docker Desktop
3. Install the Python Docker package:
   ```bash
   pip install docker
   ```

#### Linux (Ubuntu/Debian)
```bash
# Install Docker
sudo apt-get update
sudo apt-get install docker.io

# Start Docker service
sudo systemctl start docker
sudo systemctl enable docker

# Add your user to docker group (optional, to run without sudo)
sudo usermod -aG docker $USER
# Log out and back in for this to take effect

# Install Python Docker package
pip install docker
```

#### Linux (RedHat/CentOS/Fedora)
```bash
# Install Docker
sudo yum install docker

# Start Docker service
sudo systemctl start docker
sudo systemctl enable docker

# Add your user to docker group (optional)
sudo usermod -aG docker $USER
# Log out and back in for this to take effect

# Install Python Docker package
pip install docker
```

#### Windows
1. Download [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop)
2. Install and start Docker Desktop
3. Install the Python Docker package:
   ```bash
   pip install docker
   ```

### Verifying Docker Installation

After installation, verify Docker is running:
```bash
docker --version
docker ps
```

## Option 2: Local Octave Installation

If you prefer not to use Docker, you can install Octave locally.

### Installing Octave

#### macOS
```bash
# Using Homebrew (recommended)
brew install octave

# Install Python interface
pip install oct2py
```

#### Linux (Ubuntu/Debian)
```bash
# Install Octave
sudo apt-get update
sudo apt-get install octave

# Install Python interface
pip install oct2py
```

#### Linux (RedHat/CentOS/Fedora)
```bash
# Install Octave
sudo yum install octave

# Install Python interface
pip install oct2py
```

#### Windows
1. Download Octave from [GNU Octave website](https://www.gnu.org/software/octave/)
2. Install Octave and make sure it's in your PATH
3. Install Python interface:
   ```bash
   pip install oct2py
   ```

### Verifying Octave Installation

After installation, verify Octave is available:
```bash
octave --version
# or
octave-cli --version
```

## How BasisREMY Uses Octave

### Automatic Selection
When you run a simulation, BasisREMY:

1. **Checks for existing Octave instance** - If you've already run a simulation in this session, it reuses the existing instance
2. **Tries Docker first** - Looks for Docker daemon and pulls the Octave image if needed
3. **Falls back to local Octave** - If Docker is not available, tries to use local installation
4. **Shows instructions** - If neither is available, displays a helpful dialog with installation instructions

### Preference Control

You can control the preference in your backend code:
```python
# Prefer Docker (default)
backend.initialize_octave(prefer_docker=True)

# Prefer local Octave
backend.initialize_octave(prefer_docker=False)
```

## Troubleshooting

### Docker Issues

**Problem**: "Failed to connect to Docker"
- **Solution**: Make sure Docker Desktop is running
- On Linux: `sudo systemctl start docker`
- On macOS/Windows: Start Docker Desktop app

**Problem**: "Permission denied" when accessing Docker (Linux)
- **Solution**: Add your user to the docker group:
  ```bash
  sudo usermod -aG docker $USER
  ```
  Then log out and back in.

**Problem**: Docker image download fails
- **Solution**: Check your internet connection
- The Octave Docker image is downloaded automatically on first use

### Local Octave Issues

**Problem**: "octave: command not found"
- **Solution**: Make sure Octave is installed and in your PATH
- Try: `which octave` or `which octave-cli`

**Problem**: "Failed to initialize local Octave"
- **Solution**: Make sure oct2py is installed:
  ```bash
  pip install oct2py
  ```

**Problem**: oct2py import errors
- **Solution**: Some systems may need additional packages:
  ```bash
  pip install pyzmq
  ```

## For Developers

### Backend Implementation

Backends that require Octave should:

1. Set `self.requires_octave = True` in `__init__`
2. Don't initialize Octave in `__init__` (lazy initialization)
3. Call `self.initialize_octave()` in `run_simulation()` if `self.octave is None`
4. Implement `setup_octave_paths()` to configure Octave environment

Example:
```python
class MyBackend(Backend):
    def __init__(self):
        super().__init__()
        self.name = 'MyBackend'
        self.requires_octave = True  # Mark as requiring Octave
        # Don't initialize self.octave here
    
    def setup_octave_paths(self):
        """Configure Octave paths for this backend"""
        self.octave.addpath('./my/toolbox/path/')
    
    def run_simulation(self, params, progress_callback=None):
        # Initialize Octave on first use
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
            self.setup_octave_paths()
        
        # Now use self.octave...
```

### OctaveManager API

The `OctaveManager` class provides:
- `check_docker_availability()` - Returns True if Docker is available
- `check_local_octave_availability()` - Returns True if local Octave is available  
- `initialize_octave(prefer_docker=True)` - Initialize and return Octave instance
- `get_runtime_info()` - Get information about available runtimes

## Frequently Asked Questions

**Q: Which option should I use?**
A: Docker is recommended for most users as it's easier to set up and more consistent.

**Q: Can I switch between Docker and local Octave?**
A: Yes! BasisREMY will automatically use whichever is available. If both are available, it prefers Docker by default.

**Q: Do I need to download the Octave Docker image manually?**
A: No, BasisREMY automatically downloads it on first use if Docker is available.

**Q: How much disk space does the Docker option need?**
A: The Octave Docker image is approximately 1-2 GB.

**Q: Will this work offline?**
A: Yes, once the Docker image is downloaded or Octave is installed locally, you can work offline.

**Q: What if I already have Octave installed?**
A: BasisREMY will detect and use it! You can choose to use Docker instead if you prefer.
