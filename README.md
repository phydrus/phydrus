<img src=https://github.com/phydrus/phydrus/blob/master/docs/_static/logo.png width=120, align=left>

# Phydrus: Python implementation of HYDRUS-1D

<a href="http://www.gnu.org/licenses/gpl-3.0.txt"><img src=https://img.shields.io/github/license/phydrus/phydrus> </a>
<a href="https://pypi.python.org/pypi/phydrus"> <img src=https://img.shields.io/pypi/pyversions/phydrus> </a>
<a href="https://github.com/pastas/phydrus/releases"> <img src=https://img.shields.io/github/release-pre/phydrus/phydrus> </a>
<a href="https://phydrus.readthedocs.io/en/latest/?badge=latest"> <img src="https://readthedocs.org/projects/phydrus/badge/?version=latest"></a>
[![Build Status](https://travis-ci.org/phydrus/phydrus.svg?branch=master)](https://travis-ci.org/raoulcollenteur/phydrus)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4235a7486bea41c8b09e2acfa5e93e5f)](https://www.codacy.com/gh/phydrus/phydrus?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=phydrus/phydrus&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://api.codacy.com/project/badge/Coverage/4235a7486bea41c8b09e2acfa5e93e5f)](https://www.codacy.com/gh/phydrus/phydrus?utm_source=github.com&utm_medium=referral&utm_content=phydrus/phydrus&utm_campaign=Badge_Coverage)
<br>
This package provides a Python implementation of the HYDRUS-1D unsaturated zone model developed by Šimůnek, J., M. Th. van Genuchten, and M. Šejna.
More information on the HYDRUS-1D model is available [here](https://www.pc-progress.com/en/Default.aspx?hydrus-1d). This software is licenced under the GNU GENERAL PUBLIC LICENSE found [here](http://www.gnu.org/licenses/gpl-3.0.txt). The Phydrus code is developed by R.A. Collenteur and M. Vremec. With Phydrus, a HYDRUS-1D model can be created, calibrated and visualized through Python scripts, making it easy to adjust the model and providing a 100% reproducible workflow of your modeling process.

## Examples and Documentation
Examples of using Phydrus can be found in the example folder. This folder also contains a number of Jupyter Notebooks that thoroughly explain the use of the software. Documentation is hosted on ReadTheDocs: https://phydrus.readthedocs.io. 

## Bug reports and Questions
Phydrus is in active development, and bug reports are welcome as [GitHub Issues](https://github.com/phydrus/phydrus/issues). General questions or 
discussions are possible through [GitHub Discussions](https://github.com/phydrus/phydrus/discussions)

## Installing Phydrus
### 1. Installing the Python package
The Phydrus package is available from the Pypi package index and can be installed as follows:

`>>> pip install phydrus`

To install in developer mode, use the following syntax:

`>>> pip install -e .`

### 2. Compiling the source code
Before you can use Phydrus the adapted Fortran77 files need to be compiled to an executable. Phydrus-optimized Fortran Hydrus-1D code is available 
from a [dedicated repository](https://github.com/phydrus/source_code). Detailed instructions on how to compile the Fortran code can be found on the [Documentation website](https://phydrus.readthedocs.io/en/latest/getting_started.html).
 
## Developing Phydrus
Phydrus is a community effort and help is always welcome. If you have found a bug, please open a GitHub issue to report it. Pull requests including bug fixes and new functionality are very welcome and will be accepted on the Dev-branch of this repository.

## Citing Phydrus
If you use phydrus for one of your projects, we ask that you cite the code as follows:
*Collenteur, R.A., Brunetti, G., and M. Vremec (2019) Phydrus: Python implementation of the HYDRUS-1D unsaturated zone model. Version X.X.X* 
