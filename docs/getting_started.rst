Getting started with Phydrus
============================

Installing Python
-----------------
To install Phydrus, a working version of Python 3.6 or higher has to be
installed on your computer. We recommend using the `Anaconda Distribution
<https://www.continuum.io/downloads>`_ of Python. This Python distribution
includes most of the python package dependencies and the Jupyter Notebook
software to run the notebooks. Moreover, it includes the Graphical User
Interface (GUI) Spyder to start scripting in Python. However, you are free
to install any Python distribution you want.

Installing the Python package
-----------------------------
The Phydrus package will be available on the Pypi package index as the
software moves towards production ready software. To install in developer
mode, use the following syntax:

>>> pip install -e .

Compiling the source code
-------------------------
Before you can use Phydrus the adapted Fortran77 files need to be compiled to
an executable. The makefile (MacOS/Linux) and the make.bat (Windows) in the
`source` folder are available to create the Hydrus executable. Move into the
source folder (which you can download from this GitHub page) and use the
following syntax in your terminal or windows command line to compile the
source code:

>>> make

This should create a Windows or Linux Executable that can be used to run the
HYDRUS-1D simulation. In the Python code, you have to reference to the
location of the executable so you can store it anywhere you want. For MacOS,
you may need to install gfortran to compile the source code.

Dependencies
------------
Phydrus depends on a number of Python packages, of which all of the necessary
are automatically installed when using the pip install manager. To
summarize, the following packages are necessary for a minimal function
installation of Phydrus:

.. include:: ../requirements.txt
    :literal:

Updating Phydrus
----------------
If you have already installed Phydrus, it is possible to update Phydrus
easily. To update, open a Windows command screen or a Mac terminal and type::

>>> pip install phydrus --upgrade
