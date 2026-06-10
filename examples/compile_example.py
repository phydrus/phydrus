"""
Example script demonstrating the new compilation functionality in Phydrus.

This script shows how to:
1. Check for available compilers
2. Find existing executables
3. Compile HYDRUS-1D from source code
4. Use the compiled executable in a model

"""

import os
import phydrus as ps


def main():
    print("Phydrus Compilation Example")
    print("=" * 40)
    
    # 1. Check compiler information
    print("\n1. Checking available compilers...")
    compiler_info = ps.get_compiler_info()
    print(f"System: {compiler_info['system']} ({compiler_info['machine']})")
    print(f"pymake available: {compiler_info['pymake_available']}")
    print(f"gfortran available: {compiler_info['gfortran_available']}")
    print(f"make available: {compiler_info['make_available']}")
    print(f"Available compilers: {', '.join(compiler_info['compilers']) if compiler_info['compilers'] else 'None'}")
    
    # 2. Try to find existing executable
    print("\n2. Looking for existing HYDRUS-1D executable...")
    existing_exe = ps.find_executable()
    if existing_exe:
        print(f"Found existing executable: {existing_exe}")
    else:
        print("No existing executable found.")
    
    # 3. Get default executable path
    print("\n3. Default executable path...")
    default_path = ps.get_default_executable_path()
    print(f"Default path: {default_path}")
    
    # 4. Example: Compile HYDRUS-1D (commented out by default)
    print("\n4. Compilation example (commented out):")
    print("# To compile HYDRUS-1D, you can use:")
    print("# exe_path = ps.compile_hydrus()")
    print("# This will download source code and compile it automatically")
    print("#")
    print("# Or to compile to a specific location:")
    print("# exe_path = ps.compile_hydrus(target_exe='/path/to/hydrus')")
    print("#")
    print("# You can also specify the compilation method:")
    print("# exe_path = ps.compile_hydrus(method='make')  # Use make")
    print("# exe_path = ps.compile_hydrus(method='pymake')  # Use pymake")
    
    # 5. Example: Ensure executable is available
    print("\n5. Ensuring executable is available...")
    try:
        exe_path = ps.ensure_executable(compile_if_missing=False)
        print(f"Executable available at: {exe_path}")
    except FileNotFoundError:
        print("No executable found. To automatically compile, use:")
        print("exe_path = ps.ensure_executable(compile_if_missing=True)")
    
    # 6. Example: Using with Model class
    print("\n6. Using with Model class...")
    if existing_exe:
        # Use existing executable
        ws = "compile_example"
        ml = ps.Model(exe_name=existing_exe, ws_name=ws, name="compilation_example")
        print(f"Model created with executable: {ml.exe_name}")
    else:
        print("To create a model with automatic compilation if executable is missing:")
        print("ws = 'compile_example'")
        print("exe = 'hydrus'  # or path to executable")
        print("ml = ps.Model(exe_name=exe, ws_name=ws, compile_if_missing=True)")
    
    print("\nExample completed!")


if __name__ == "__main__":
    main()