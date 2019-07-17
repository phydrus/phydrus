# Pydrus
<a href="http://www.gnu.org/licenses/gpl-3.0.txt"><img src=https://img.shields.io/pypi/l/pastas.svg></a>

This package provides a Python implementation of the HYDRUS-1D unsaturated zone model. More information on the HYDRUS-1D model is available [here](https://www.pc-progress.com/en/Default.aspx?hydrus-1d). This software is licenced under the GNU GENERAL PUBLIC LICENSE found [here](http://www.gnu.org/licenses/gpl-3.0.txt). Pydrus 

# Examples and Documentation
Examples of using Pydrus can be found in the example folder. This folder also contains a number of Jupyter Notebooks that thoroughly explain the use of the software. Documentation is still to be developed.

# Installing the Pydrus Python package
The Pydrus package will be available on the Pypi package index as the software moves towards production ready software. To install in developer mode, use the following syntax:

`>>> pip install -e .`

# Compiling the source code
Before you can use Pydrus the adapted Fortran77 files need to be compiled to an executable. The makefile (MacOS/Linux) and the make.bat (Windows) in the `source` folder are available to create the hydrus executable. Move into the source fodler and use the following syntax in your terminal or windows command line to compile the source code:
 
`>>> make`
 
This should create a Windows or Linux Executable that can be used to run the HYDRUS-1D simulation.  
 
# Developing Pydrus
Pydrus is a community effort and help is always welcome. If you have found abug, please open a GitHub issue to report it. Pull requests including bug fixes and new functionality are veery welcome and will be accepted on the Dev-branch of this repository.
