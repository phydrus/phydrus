"""
Tests for the compilation module.

These tests verify the functionality of the phydrus.compile module.
"""

import os
import tempfile
import shutil
import pytest
import platform
from unittest.mock import patch, MagicMock

import phydrus as ps
from phydrus.compile import (
    get_default_executable_path, find_executable, get_compiler_info,
    CompilationError, DownloadError
)


class TestExecutablePaths:
    """Test functions related to finding and determining executable paths."""
    
    def test_get_default_executable_path(self):
        """Test getting default executable path."""
        path = get_default_executable_path()
        assert isinstance(path, str)
        assert len(path) > 0
        
        # Check that it ends with appropriate executable name for the platform
        if platform.system() == "Windows":
            assert path.endswith("hydrus.exe") or path.endswith("hydrus1d.exe")
        else:
            assert path.endswith("hydrus") or path.endswith("hydrus1d")
    
    def test_find_executable_with_existing_exe(self):
        """Test finding executable when one exists in current directory."""
        # Create a temporary directory with a fake executable
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            os.chdir(tmpdir)
            
            try:
                # Create a fake executable
                if platform.system() == "Windows":
                    exe_name = "hydrus.exe"
                else:
                    exe_name = "hydrus"
                
                with open(exe_name, 'w') as f:
                    f.write("fake executable")
                
                # Make it executable on Unix-like systems
                if platform.system() != "Windows":
                    os.chmod(exe_name, 0o755)
                
                # Test finding it
                found_path = find_executable()
                assert found_path is not None
                assert os.path.exists(found_path)
                
            finally:
                os.chdir(original_cwd)
    
    def test_find_executable_nonexistent(self):
        """Test finding executable when none exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            original_path = os.environ.get("PATH", "")
            os.chdir(tmpdir)
            
            try:
                # Clear PATH to ensure no existing executables are found
                os.environ["PATH"] = tmpdir
                
                # Test that no executable is found
                found_path = find_executable()
                assert found_path is None
                
            finally:
                os.chdir(original_cwd)
                os.environ["PATH"] = original_path


class TestCompilerInfo:
    """Test compiler information functions."""
    
    def test_get_compiler_info_structure(self):
        """Test that get_compiler_info returns expected structure."""
        info = get_compiler_info()
        
        assert isinstance(info, dict)
        assert "system" in info
        assert "machine" in info
        assert "pymake_available" in info
        assert "gfortran_available" in info
        assert "make_available" in info
        assert "compilers" in info
        
        # Check that system and machine are strings
        assert isinstance(info["system"], str)
        assert isinstance(info["machine"], str)
        
        # Check that boolean fields are booleans
        assert isinstance(info["pymake_available"], bool)
        assert isinstance(info["gfortran_available"], bool)
        assert isinstance(info["make_available"], bool)
        
        # Check that compilers is a list
        assert isinstance(info["compilers"], list)


class TestCompilationFunctions:
    """Test compilation functions with mocking to avoid actual compilation."""
    
    @patch('phydrus.compile.download_source_code')
    @patch('phydrus.compile.compile_with_make')
    def test_compile_hydrus_with_make(self, mock_compile_make, mock_download):
        """Test compile_hydrus function using make method."""
        # Setup mocks
        mock_download.return_value = "/fake/source/dir"
        mock_compile_make.return_value = "/fake/executable"
        
        # Call the function
        result = ps.compile_hydrus(method="make")
        
        # Verify it was called with correct arguments
        assert mock_download.called
        assert mock_compile_make.called
        assert result == "/fake/executable"
    
    @patch('phydrus.compile.download_source_code')
    @patch('phydrus.compile.compile_with_pymake')
    def test_compile_hydrus_with_pymake(self, mock_compile_pymake, mock_download):
        """Test compile_hydrus function using pymake method."""
        # Setup mocks
        mock_download.return_value = "/fake/source/dir"
        mock_compile_pymake.return_value = "/fake/executable"
        
        # Call the function
        result = ps.compile_hydrus(method="pymake")
        
        # Verify it was called with correct arguments
        assert mock_download.called
        assert mock_compile_pymake.called
        assert result == "/fake/executable"
    
    @patch('phydrus.compile.download_source_code')
    @patch('phydrus.compile.compile_with_pymake')
    @patch('phydrus.compile.compile_with_make')
    def test_compile_hydrus_auto_fallback(self, mock_compile_make, mock_compile_pymake, mock_download):
        """Test compile_hydrus function with auto method falling back from pymake to make."""
        # Setup mocks
        mock_download.return_value = "/fake/source/dir"
        mock_compile_pymake.side_effect = CompilationError("pymake failed")
        mock_compile_make.return_value = "/fake/executable"
        
        # Call the function
        result = ps.compile_hydrus(method="auto")
        
        # Verify fallback behavior
        assert mock_download.called
        assert mock_compile_pymake.called
        assert mock_compile_make.called
        assert result == "/fake/executable"


class TestEnsureExecutable:
    """Test ensure_executable function."""
    
    def test_ensure_executable_with_existing(self):
        """Test ensure_executable when executable already exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a fake executable
            if platform.system() == "Windows":
                exe_path = os.path.join(tmpdir, "hydrus.exe")
            else:
                exe_path = os.path.join(tmpdir, "hydrus")
            
            with open(exe_path, 'w') as f:
                f.write("fake executable")
            
            if platform.system() != "Windows":
                os.chmod(exe_path, 0o755)
            
            # Test that existing executable is returned
            result = ps.ensure_executable(exe_path, compile_if_missing=False)
            assert result == exe_path
    
    @patch('phydrus.compile.compile_hydrus')
    def test_ensure_executable_compile_missing(self, mock_compile):
        """Test ensure_executable when executable is missing and compilation is requested."""
        mock_compile.return_value = "/fake/compiled/executable"
        
        # Test with non-existent path
        result = ps.ensure_executable("/nonexistent/path", compile_if_missing=True)
        
        # Verify compilation was triggered
        assert mock_compile.called
        assert result == "/fake/compiled/executable"
    
    def test_ensure_executable_no_compile_missing(self):
        """Test ensure_executable when executable is missing and compilation is not requested."""
        with pytest.raises(FileNotFoundError):
            ps.ensure_executable("/nonexistent/path", compile_if_missing=False)


class TestModelIntegration:
    """Test integration of compilation functionality with Model class."""
    
    def test_model_set_executable_with_existing(self):
        """Test Model.set_executable with existing executable."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a fake executable
            if platform.system() == "Windows":
                exe_path = os.path.join(tmpdir, "hydrus.exe")
            else:
                exe_path = os.path.join(tmpdir, "hydrus")
            
            with open(exe_path, 'w') as f:
                f.write("fake executable")
            
            if platform.system() != "Windows":
                os.chmod(exe_path, 0o755)
            
            # Create model and set executable
            ml = ps.Model(exe_name=exe_path, ws_name=os.path.join(tmpdir, "workspace"))
            assert ml.exe_name == exe_path
    
    @patch('phydrus.compile.ensure_executable')
    def test_model_set_executable_compile_if_missing(self, mock_ensure):
        """Test Model.set_executable with compile_if_missing=True."""
        mock_ensure.return_value = "/fake/compiled/executable"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            exe_path = os.path.join(tmpdir, "nonexistent.exe")
            
            # Create model with compile_if_missing
            ml = ps.Model(exe_name=exe_path, ws_name=os.path.join(tmpdir, "workspace"), 
                         compile_if_missing=True)
            
            # Verify ensure_executable was called
            assert mock_ensure.called


class TestCommandLineInterface:
    """Test command line interface functionality."""
    
    @patch('phydrus.compile.compile_hydrus')
    @patch('sys.argv', ['phydrus-compile', '--info'])
    def test_main_info(self, mock_compile):
        """Test main function with --info flag."""
        from phydrus.compile import main
        
        # This should not raise an exception
        main()
        
        # compile_hydrus should not be called
        assert not mock_compile.called
    
    @patch('phydrus.compile.compile_hydrus')
    @patch('sys.argv', ['phydrus-compile'])
    def test_main_default(self, mock_compile):
        """Test main function with default arguments."""
        from phydrus.compile import main
        
        mock_compile.return_value = "/fake/executable"
        
        # This should not raise an exception
        main()
        
        # compile_hydrus should be called
        assert mock_compile.called


# Test data for mocking
@pytest.fixture
def mock_source_dir():
    """Create a mock source directory for testing."""
    tmpdir = tempfile.mkdtemp()
    
    # Create a simple Fortran file
    with open(os.path.join(tmpdir, "test.f"), 'w') as f:
        f.write("      PROGRAM TEST\n      END\n")
    
    yield tmpdir
    
    # Cleanup
    shutil.rmtree(tmpdir, ignore_errors=True)


class TestDownloadSourceCode:
    """Test source code download functionality."""
    
    @patch('urllib.request.urlretrieve')
    @patch('zipfile.ZipFile')
    def test_download_source_code(self, mock_zipfile, mock_urlretrieve):
        """Test download_source_code function."""
        # Setup mocks
        mock_urlretrieve.return_value = None
        mock_zip = MagicMock()
        mock_zipfile.return_value.__enter__.return_value = mock_zip
        
        # Mock the zip file contents
        mock_zip.namelist.return_value = ['source_code-master/', 'source_code-master/test.f']
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # This should not raise an exception with our mocks
            result = ps.download_source_code(target_dir=tmpdir)
            
            # Verify the download was attempted
            assert mock_urlretrieve.called
            assert mock_zipfile.called