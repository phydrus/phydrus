"""
Compilation utilities for Phydrus HYDRUS-1D executable.

This module provides functionality to download and compile the HYDRUS-1D
Fortran source code from the phydrus/source_code repository.

"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from logging import getLogger

logger = getLogger(__name__)


class CompilationError(Exception):
    """Exception raised when compilation fails."""

    pass


class DownloadError(Exception):
    """Exception raised when downloading source code fails."""

    pass


def get_default_executable_path():
    """
    Get the default path where the HYDRUS-1D executable should be located.

    Returns
    -------
    str
        Default path for the executable.
    """
    # Check if there's an executable in the current directory
    if platform.system() == "Windows":
        default_names = ["hydrus.exe", "hydrus1d.exe"]
    else:
        default_names = ["hydrus", "hydrus1d"]

    for name in default_names:
        if os.path.exists(name):
            return os.path.abspath(name)

    # Check in the examples directory
    examples_dir = os.path.join(os.path.dirname(__file__), "..", "examples")
    for name in default_names:
        exe_path = os.path.join(examples_dir, name)
        if os.path.exists(exe_path):
            return os.path.abspath(exe_path)

    # Return a default path in the user's home directory
    home_dir = os.path.expanduser("~")
    if platform.system() == "Windows":
        return os.path.join(home_dir, "hydrus", "hydrus.exe")
    else:
        return os.path.join(home_dir, ".local", "bin", "hydrus")


def find_executable():
    """
    Try to find the HYDRUS-1D executable in common locations.

    Returns
    -------
    str or None
        Path to the executable if found, None otherwise.
    """
    # Check PATH environment variable
    if platform.system() == "Windows":
        exe_names = ["hydrus.exe", "hydrus1d.exe"]
    else:
        exe_names = ["hydrus", "hydrus1d"]

    for path in os.environ.get("PATH", "").split(os.pathsep):
        for exe_name in exe_names:
            full_path = os.path.join(path, exe_name)
            if os.path.exists(full_path) and os.access(full_path, os.X_OK):
                return full_path

    # Check current directory and parent directories
    current_dir = os.getcwd()
    for _ in range(3):  # Check up to 3 levels up
        for exe_name in exe_names:
            full_path = os.path.join(current_dir, exe_name)
            if os.path.exists(full_path) and os.access(full_path, os.X_OK):
                return full_path
        current_dir = os.path.dirname(current_dir)

    return None


def download_source_code(target_dir=None, repo_url=None, branch="main"):
    """
    Download the phydrus/source_code repository.

    Parameters
    ----------
    target_dir : str, optional
        Directory to download the source code to. If None, creates a
        temporary directory.
    repo_url : str, optional
        URL of the source code repository. Defaults to phydrus/source_code.
    branch : str, optional
        Branch to download. Defaults to "master".

    Returns
    -------
    str
        Path to the downloaded source code directory.

    Raises
    ------
    DownloadError
        If downloading the source code fails.
    """
    if repo_url is None:
        repo_url = "https://github.com/phydrus/source_code"

    if target_dir is None:
        target_dir = tempfile.mkdtemp(prefix="phydrus_source_")
    else:
        target_dir = os.path.abspath(target_dir)

    # Create target directory if it doesn't exist
    os.makedirs(target_dir, exist_ok=True)

    # Download the repository as a zip file
    zip_url = f"{repo_url}/archive/refs/heads/{branch}.zip"
    zip_path = os.path.join(target_dir, "source_code.zip")

    logger.info(f"Downloading source code from {zip_url}")

    try:
        # Download the zip file
        urllib.request.urlretrieve(zip_url, zip_path)

        # Extract the zip file
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(target_dir)

        # Remove the zip file
        os.remove(zip_path)

        # Find the extracted directory
        extracted_dir = None
        for item in os.listdir(target_dir):
            item_path = os.path.join(target_dir, item)
            if os.path.isdir(item_path) and item.startswith("source_code"):
                extracted_dir = item_path
                break

        if extracted_dir is None:
            raise DownloadError("Could not find extracted source code directory")

        logger.info(f"Source code downloaded to {extracted_dir}")
        return extracted_dir

    except Exception as e:
        # Clean up on error
        if os.path.exists(zip_path):
            os.remove(zip_path)
        if os.path.exists(target_dir):
            shutil.rmtree(target_dir)
        raise DownloadError(f"Failed to download source code: {e}")


def compile_with_make(source_dir, target_exe=None, make_command="make"):
    """
    Compile the HYDRUS-1D source code using make.

    Parameters
    ----------
    source_dir : str
        Directory containing the source code.
    target_exe : str, optional
        Path for the output executable. If None, uses the source directory.
    make_command : str, optional
        Make command to use. Defaults to "make".

    Returns
    -------
    str
        Path to the compiled executable.

    Raises
    ------
    CompilationError
        If compilation fails.
    """
    # Store original target if provided
    original_target_exe = target_exe

    if target_exe is None:
        if platform.system() == "Windows":
            target_exe = os.path.join(source_dir, "hydrus.exe")
        else:
            target_exe = os.path.join(source_dir, "hydrus")
    else:
        target_exe = os.path.abspath(target_exe)

    # Change to source directory
    original_dir = os.getcwd()
    os.chdir(source_dir)

    try:
        logger.info(f"Compiling HYDRUS-1D in {source_dir}")
        logger.info(f"Using make command: {make_command}")

        # Run make
        result = subprocess.run(
            [make_command], capture_output=True, text=True, check=False
        )

        if result.returncode != 0:
            logger.error(f"Compilation failed with return code {result.returncode}")
            logger.error(f"STDOUT: {result.stdout}")
            logger.error(f"STDERR: {result.stderr}")
            raise CompilationError(f"Compilation failed: {result.stderr}")
        print(result.returncode)
        # Check if executable was created at the target location
        if not os.path.exists(target_exe):
            # Try alternative names
            if platform.system() == "Windows":
                alt_names = ["hydrus.exe", "hydrus1d.exe"]
            else:
                alt_names = ["hydrus", "hydrus1d"]

            found_exe = None
            for name in alt_names:
                alt_path = os.path.join(source_dir, name)
                if os.path.exists(alt_path):
                    found_exe = alt_path
                    break

            if found_exe is None:
                raise CompilationError(f"Executable not found in {source_dir}")

            # If a custom target was provided and it differs from where we found it,
            # move the executable to the desired location
            if (
                original_target_exe is not None
                and os.path.abspath(original_target_exe) != found_exe
            ):
                # Ensure target directory exists
                target_dir = os.path.dirname(target_exe)
                if target_dir and not os.path.exists(target_dir):
                    os.makedirs(target_dir, exist_ok=True)
                shutil.move(found_exe, target_exe)
                logger.info(f"Moved executable from {found_exe} to {target_exe}")
            else:
                target_exe = found_exe

        logger.info(f"Successfully compiled executable: {target_exe}")
        return target_exe

    finally:
        os.chdir(original_dir)


def compile_with_pymake(
    source_dir, target_exe=None, fc="gfortran", include_subdirs=True, subdirs=None
):
    """
    Compile the HYDRUS-1D source code using pymake.

    Parameters
    ----------
    source_dir : str
        Directory containing the source code.
    target_exe : str, optional
        Path for the output executable. If None, uses the source directory.
    fc : str, optional
        Fortran compiler to use. Defaults to "gfortran".
    include_subdirs : bool, optional
        Whether to include subdirectories. Defaults to True.
    subdirs : list, optional
        List of subdirectories to include. If None and include_subdirs is True,
        all subdirectories will be included.

    Returns
    -------
    str
        Path to the compiled executable.

    Raises
    ------
    CompilationError
        If compilation fails or pymake is not available.
    """
    try:
        import pymake
    except ImportError:
        raise CompilationError(
            "pymake is not installed. Install it with: pip install mfpymake"
        )

    if target_exe is None:
        if platform.system() == "Windows":
            target_exe = os.path.join(source_dir, "hydrus.exe")
        else:
            target_exe = os.path.join(source_dir, "hydrus")
    else:
        target_exe = os.path.abspath(target_exe)

    # Create pymake object
    pm = pymake.Pymake()
    pm.srcdir = source_dir
    pm.target = target_exe
    pm.fc = fc
    pm.include_subdirs = include_subdirs

    if subdirs is not None:
        pm.subdirs = subdirs

    logger.info(f"Compiling HYDRUS-1D with pymake using {fc} compiler")
    logger.info(f"Source directory: {source_dir}")
    logger.info(f"Target executable: {target_exe}")

    try:
        # Build the executable
        pm.build()

        # Check if executable was created
        if not os.path.exists(target_exe):
            raise CompilationError(f"Executable not found at {target_exe}")

        logger.info(f"Successfully compiled executable with pymake: {target_exe}")
        return target_exe

    except Exception as e:
        raise CompilationError(f"pymake compilation failed: {e}")


def compile_hydrus(
    source_dir=None, target_exe=None, method="auto", keep_source=False, **kwargs
):
    """
    Compile the HYDRUS-1D executable from source code.

    This is the main function for compiling the HYDRUS-1D executable.
    It can automatically download the source code from the phydrus/source_code
    repository and compile it using either make or pymake.

    Parameters
    ----------
    source_dir : str, optional
        Directory containing the source code. If None, the source code
        will be downloaded from the repository.
    target_exe : str, optional
        Path for the output executable. If None, uses a default location.
    method : str, optional
        Compilation method: "auto", "make", or "pymake".
        "auto" will try pymake first, then fall back to make.
    keep_source : bool, optional
        Whether to keep the downloaded source code after compilation.
        Defaults to False.
    **kwargs
        Additional arguments passed to the compilation functions.

    Returns
    -------
    str
        Path to the compiled executable.

    Raises
    ------
    CompilationError
        If compilation fails.
    DownloadError
        If downloading source code fails.

    Examples
    --------
    >>> from phydrus.compile import compile_hydrus
    >>> exe_path = compile_hydrus()
    >>> print(f"Executable compiled at: {exe_path}")

    >>> # Compile to a specific location
    >>> exe_path = compile_hydrus(target_exe="/path/to/hydrus")

    >>> # Use make instead of pymake
    >>> exe_path = compile_hydrus(method="make")
    """
    # Determine source directory
    temp_source_dir = None
    if source_dir is None:
        temp_source_dir = tempfile.mkdtemp(prefix="phydrus_compile_")
        source_dir = download_source_code(target_dir=temp_source_dir)
    else:
        source_dir = os.path.abspath(source_dir)

    # Determine target executable
    if target_exe is None:
        target_exe = get_default_executable_path()
    else:
        target_exe = os.path.abspath(target_exe)

    # Create target directory if it doesn't exist
    target_dir = os.path.dirname(target_exe)
    if target_dir and not os.path.exists(target_dir):
        os.makedirs(target_dir, exist_ok=True)

    try:
        if method == "auto":
            # Try pymake first, then fall back to make
            try:
                exe_path = compile_with_pymake(source_dir, target_exe, **kwargs)
            except (CompilationError, ImportError) as e:
                logger.info(f"pymake method failed, trying make: {e}")
                exe_path = compile_with_make(source_dir, target_exe, **kwargs)
        elif method == "pymake":
            exe_path = compile_with_pymake(source_dir, target_exe, **kwargs)
        elif method == "make":
            exe_path = compile_with_make(source_dir, target_exe, **kwargs)
        else:
            raise CompilationError(f"Unknown compilation method: {method}")

        return exe_path

    finally:
        # Clean up temporary source directory if requested
        if temp_source_dir and not keep_source:
            shutil.rmtree(temp_source_dir, ignore_errors=True)


def ensure_executable(exe_path=None, compile_if_missing=True, **kwargs):
    """
    Ensure that a HYDRUS-1D executable is available.

    This function checks if an executable exists at the specified path,
    and optionally compiles it if it's missing.

    Parameters
    ----------
    exe_path : str, optional
        Path to the executable. If None, tries to find an existing executable
        or uses a default location.
    compile_if_missing : bool, optional
        Whether to compile the executable if it's not found. Defaults to True.
    **kwargs
        Additional arguments passed to compile_hydrus if compilation is needed.

    Returns
    -------
    str
        Path to the executable.

    Raises
    ------
    CompilationError
        If compilation is requested but fails.
    FileNotFoundError
        If no executable is found and compilation is not requested.

    Examples
    --------
    >>> from phydrus.compile import ensure_executable
    >>> exe_path = ensure_executable()
    >>> print(f"Using executable: {exe_path}")
    """
    if exe_path is None:
        # Try to find existing executable
        exe_path = find_executable()
        if exe_path is None:
            exe_path = get_default_executable_path()
    else:
        exe_path = os.path.abspath(exe_path)

    # Check if executable exists and is executable
    if os.path.exists(exe_path) and os.access(exe_path, os.X_OK):
        logger.info(f"Found existing executable: {exe_path}")
        return exe_path

    if not compile_if_missing:
        raise FileNotFoundError(f"Executable not found at {exe_path}")

    logger.info(f"Executable not found at {exe_path}, compiling...")

    # Compile the executable
    return compile_hydrus(target_exe=exe_path, **kwargs)


def get_compiler_info():
    """
    Get information about available compilers and compilation options.

    Returns
    -------
    dict
        Dictionary with compiler information.
    """
    info = {
        "system": platform.system(),
        "machine": platform.machine(),
        "pymake_available": False,
        "gfortran_available": False,
        "make_available": False,
        "compilers": [],
    }

    # Check for pymake
    try:
        import pymake

        info["pymake_available"] = True
        info["pymake_version"] = getattr(pymake, "__version__", "unknown")
    except ImportError:
        pass

    # Check for gfortran
    try:
        result = subprocess.run(
            ["gfortran", "--version"], capture_output=True, check=False
        )
        if result.returncode == 0:
            info["gfortran_available"] = True
            info["compilers"].append("gfortran")
    except FileNotFoundError:
        pass

    # Check for make
    try:
        result = subprocess.run(["make", "--version"], capture_output=True, check=False)
        if result.returncode == 0:
            info["make_available"] = True
    except FileNotFoundError:
        pass

    # Check for other Fortran compilers
    for compiler in ["ifort", "ifx", "pgfortran", "flang"]:
        try:
            result = subprocess.run(
                [compiler, "--version"], capture_output=True, check=False
            )
            if result.returncode == 0:
                info["compilers"].append(compiler)
        except FileNotFoundError:
            pass

    return info


def main():
    """
    Command line interface for compiling HYDRUS-1D executable.

    This function is called when using the 'phydrus-compile' command.
    """
    parser = argparse.ArgumentParser(
        description="Compile HYDRUS-1D executable for Phydrus",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  phydrus-compile                    # Compile to default location
  phydrus-compile --target /path/to/hydrus  # Compile to specific location
  phydrus-compile --method make      # Use make instead of pymake
  phydrus-compile --keep-source      # Keep downloaded source code
  phydrus-compile --info             # Show compiler information
        """,
    )

    parser.add_argument(
        "--target",
        "-t",
        help="Target path for the compiled executable. Defaults to current directory.",
    )

    parser.add_argument(
        "--source-dir",
        "-s",
        help="Directory containing source code. If not provided, will download from repository.",
    )

    parser.add_argument(
        "--method",
        "-m",
        choices=["auto", "make", "pymake"],
        default="auto",
        help="Compilation method: auto (try pymake then make), make, or pymake. Default: auto",
    )

    parser.add_argument(
        "--fc",
        "--fortran-compiler",
        default="gfortran",
        help="Fortran compiler to use (for pymake method). Default: gfortran",
    )

    parser.add_argument(
        "--keep-source",
        "-k",
        action="store_true",
        help="Keep the downloaded source code after compilation.",
    )

    parser.add_argument(
        "--info", "-i", action="store_true", help="Show compiler information and exit."
    )

    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show verbose output."
    )

    args = parser.parse_args()

    # Set up logging
    if args.verbose:
        import logging

        logging.basicConfig(level=logging.INFO)
    else:
        import logging

        logging.basicConfig(level=logging.WARNING)

    if args.info:
        info = get_compiler_info()
        print("Compiler Information:")
        print(f"  System: {info['system']} ({info['machine']})")
        print(f"  pymake available: {info['pymake_available']}")
        if info.get("pymake_version"):
            print(f"  pymake version: {info['pymake_version']}")
        print(f"  gfortran available: {info['gfortran_available']}")
        print(f"  make available: {info['make_available']}")
        print(
            f"  Available compilers: {', '.join(info['compilers']) if info['compilers'] else 'None'}"
        )
        return

    try:
        # Compile the executable
        exe_path = compile_hydrus(
            source_dir=args.source_dir,
            target_exe=args.target,
            method=args.method,
            keep_source=args.keep_source,
            fc=args.fc,
        )

        print(f"Successfully compiled HYDRUS-1D executable: {exe_path}")

    except (CompilationError, DownloadError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
