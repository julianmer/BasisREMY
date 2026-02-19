"""
Tests for core.octave_manager module
"""

import pytest
import os
from core.octave_manager import OctaveManager


@pytest.mark.core
@pytest.mark.unit
class TestOctaveManager:
    """Test OctaveManager functionality"""

    def test_init_default(self):
        """Test OctaveManager initialization with defaults"""
        # Ensure env variable is not set
        if 'BASISREMY_VERBOSE' in os.environ:
            del os.environ['BASISREMY_VERBOSE']

        manager = OctaveManager()
        assert manager.runtime_type is None
        assert manager.octave_instance is None
        assert manager.verbose == False

    def test_init_verbose(self):
        """Test OctaveManager initialization with verbose"""
        manager = OctaveManager(verbose=True)
        assert manager.verbose == True

    def test_init_verbose_from_env(self):
        """Test OctaveManager initialization with env variable"""
        os.environ['BASISREMY_VERBOSE'] = '1'
        manager = OctaveManager()
        assert manager.verbose == True
        # Clean up
        del os.environ['BASISREMY_VERBOSE']
        assert manager.verbose == True
        os.environ.pop('BASISREMY_VERBOSE', None)

    def test_check_docker_availability(self):
        """Test Docker availability check"""
        manager = OctaveManager()
        result = manager.check_docker_availability()
        assert isinstance(result, bool)

    def test_check_local_octave_availability(self):
        """Test local Octave availability check"""
        manager = OctaveManager()
        result = manager.check_local_octave_availability()
        assert isinstance(result, bool)

    @pytest.mark.requires_docker
    @pytest.mark.docker
    def test_initialize_docker(self, cleanup_docker_processes):
        """Test Docker Octave initialization"""
        manager = OctaveManager(verbose=True)
        octave = manager.initialize_octave(prefer_docker=True)

        assert octave is not None
        assert manager.runtime_type == 'docker'
        assert manager.octave_instance is not None

    @pytest.mark.requires_octave
    def test_initialize_local(self):
        """Test local Octave initialization"""
        manager = OctaveManager(verbose=False)
        octave = manager.initialize_octave(prefer_docker=False)

        assert octave is not None
        assert manager.runtime_type == 'local'
        assert manager.octave_instance is not None

    def test_get_runtime_info(self):
        """Test get_runtime_info method"""
        manager = OctaveManager(verbose=True)
        info = manager.get_runtime_info()

        assert isinstance(info, dict)
        assert 'runtime_type' in info
        assert 'docker_available' in info
        assert 'local_octave_available' in info
        assert 'verbose' in info
        assert info['verbose'] == True

    def test_initialize_with_neither_available(self):
        """Test initialization fails gracefully when nothing available"""
        # This test only makes sense if neither Docker nor Octave is available
        # Skip it if either is available
        manager = OctaveManager()
        if manager.check_docker_availability() or manager.check_local_octave_availability():
            pytest.skip("Docker or Octave is available - cannot test failure case")

        with pytest.raises(RuntimeError) as exc_info:
            manager.initialize_octave()

        assert "Octave Runtime Not Available" in str(exc_info.value)


