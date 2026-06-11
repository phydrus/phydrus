from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

# Get the version.
version = {}
with open("phydrus/version.py") as fp:
    exec(fp.read(), version)

setup(
    name='phydrus',
    version=version['__version__'],
    description='Python implementation of the HYDRUS-1D model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/phydrus/phydrus',
    license='GNU General public version 3.0',
    author='Raoul Collenteur, Matevz Vremec, Giuseppe Brunetti',
    author_email='raoul.collenteur@uni-graz.at',
    project_urls={
        'Source': 'https://github.com/phydrus/phydrus',
        'Tracker': 'https://github.com/phydrus/phydrus/issues',
        'Help': 'https://github.com/phydrus/phydrus/discussions'
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Other Audience',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering :: Hydrology',
    ],
    platforms='Windows, Mac OS-X, Linux',
    install_requires=['numpy>=1.21', 'matplotlib>=3.1', 'pandas>=2.2.0,<4.0',
                      'scipy>=1.7'],
    packages=find_packages(exclude=[]),
    package_data={"source": ["hydrus", "hydrus.exe"], },
)
