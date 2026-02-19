"""
Pytest fixtures and configuration for BasisREMY tests
"""

import os
import sys
import pytest

# Add project root to path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)


@pytest.fixture(scope="session")
def project_root_dir():
    """Return the project root directory"""
    return project_root


@pytest.fixture(scope="session")
def example_data_dir(project_root_dir):
    """Return the example data directory"""
    return os.path.join(project_root_dir, 'example_data')


@pytest.fixture(scope="session")
def output_dir(project_root_dir):
    """Return the output directory"""
    output_path = os.path.join(project_root_dir, 'output')
    os.makedirs(output_path, exist_ok=True)
    return output_path


@pytest.fixture(scope="function")
def test_output_dir(output_dir, request):
    """Create a test-specific output directory"""
    test_name = request.node.name
    test_output = os.path.join(output_dir, 'test_results', test_name)
    os.makedirs(test_output, exist_ok=True)
    return test_output


@pytest.fixture(scope="session")
def docker_available():
    """Check if Docker is available"""
    try:
        import docker
        client = docker.from_env()
        client.ping()
        return True
    except:
        return False


@pytest.fixture(scope="session")
def octave_available():
    """Check if local Octave is available"""
    import shutil
    octave_cmd = shutil.which('octave-cli') or shutil.which('octave')
    return octave_cmd is not None


@pytest.fixture(scope="function")
def enable_verbose():
    """Enable verbose mode for debugging"""
    original = os.environ.get('BASISREMY_VERBOSE', '')
    os.environ['BASISREMY_VERBOSE'] = '1'
    yield
    if original:
        os.environ['BASISREMY_VERBOSE'] = original
    else:
        os.environ.pop('BASISREMY_VERBOSE', None)


@pytest.fixture(scope="function")
def cleanup_docker_processes():
    """Cleanup Docker Octave processes after test"""
    yield
    try:
        import docker
        client = docker.from_env()
        try:
            container = client.containers.get('octave_runner')
            container.exec_run("pkill -9 octave-cli")
        except:
            pass
    except:
        pass


# Test data fixtures
@pytest.fixture(scope="session")
def philips_spar_file(example_data_dir):
    """Philips SPAR test file"""
    return os.path.join(example_data_dir, 'BigGABA_P1P_S01', 'S01_PRESS_35_act.SPAR')


@pytest.fixture(scope="session")
def siemens_rda_file(example_data_dir):
    """Siemens RDA test file (if available)"""
    # Look for .rda files in example_data
    for root, dirs, files in os.walk(example_data_dir):
        for file in files:
            if file.lower().endswith('.rda'):
                return os.path.join(root, file)
    return None


@pytest.fixture(scope="session")
def ge_p_file(example_data_dir):
    """GE P-file test file"""
    ge_dir = os.path.join(example_data_dir, 'BigGABA_G1P_S01')
    if os.path.exists(ge_dir):
        for file in os.listdir(ge_dir):
            if file.endswith('.7'):
                return os.path.join(ge_dir, file)
    return None


@pytest.fixture(scope="session")
def bruker_dat_file(example_data_dir):
    """Bruker .dat test file"""
    bruker_dir = os.path.join(example_data_dir, 'BigGABA_S1P_S01')
    if os.path.exists(bruker_dir):
        for file in os.listdir(bruker_dir):
            if file.endswith('.dat'):
                return os.path.join(bruker_dir, file)
    return None


@pytest.fixture(scope="session")
def pulse_file(project_root_dir):
    """Path to the pulse file for sLaser simulations"""
    return os.path.join(project_root_dir, 'externals', 'jbss', 'my_pulse', 'standardized_goia.txt')


# Skip markers for missing dependencies
def pytest_configure(config):
    """Configure pytest with custom markers"""
    config.addinivalue_line(
        "markers", "requires_docker: mark test as requiring Docker"
    )
    config.addinivalue_line(
        "markers", "requires_octave: mark test as requiring local Octave"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add skip markers based on availability"""
    # Check Docker availability (supports both Docker Desktop and OrbStack on macOS)
    try:
        import docker
        # Try default connection first
        try:
            client = docker.from_env()
            client.ping()
            docker_available = True
        except:
            # On macOS with OrbStack, try OrbStack socket
            import os
            orbstack_socket = os.path.expanduser('~/.orbstack/run/docker.sock')
            if os.path.exists(orbstack_socket):
                try:
                    client = docker.DockerClient(base_url=f'unix://{orbstack_socket}')
                    client.ping()
                    docker_available = True
                except:
                    docker_available = False
            else:
                docker_available = False
    except:
        docker_available = False

    import shutil
    octave_cmd = shutil.which('octave-cli') or shutil.which('octave')
    octave_available = octave_cmd is not None

    for item in items:
        if "requires_docker" in item.keywords and not docker_available:
            item.add_marker(pytest.mark.skip(reason="Docker not available"))
        if "requires_octave" in item.keywords and not octave_available:
            item.add_marker(pytest.mark.skip(reason="Local Octave not available"))

