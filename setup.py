from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

# Get the version.
version = {}
with open("pydrus/version.py") as fp:
    exec(fp.read(), version)

setup(
    name='pydrus',
    version=version['__version__'],
    description='Python implementation of the HYDRUS-1D model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/raoulcollenteur/pydrus',
    license='GNU General public version 3.0',
    author='Raoul Collenteur, Giuseppe Brunetti',
    author_email='raoul.collenteur@uni-graz.at, giuseppe.brunetti@boku.ac.at',
    project_urls={
        'Source': 'https://github.com/raoulcollenteur/pydrus',
        'Tracker': 'https://github.com/raoulcollenteur/pydrus/issues',
        'Help': 'https://stackoverflow.com/questions/tagged/pydrus'
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Other Audience',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Hydrology',
    ],
    platforms='Windows, Mac OS-X',
    install_requires=['numpy>=1.15', 'matplotlib>=2.0', 'pandas>=0.23',
                      'scipy>=1.1'],
    packages=find_packages(exclude=[]),
    package_data={"source": ["hydrus", "hydrus.exe"], },
)