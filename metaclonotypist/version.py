from setuptools import find_packages

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ""  # use '' for first of series, number for 1 and above
# _version_extra = 'dev'
_version_extra = ""  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering",
]

# Description should be a one-liner:
description = "Metaclonotype discovery pipeline"
# Long description will go up on the pypi page
long_description = """
Metaclonotypist is a fast and scalable pipeline for metaclonotype discovery build on top of pyrepseq.

Metaclonotypist can be installed from its [Github](https://github.com/qimmuno/metaclonotypist) source, by running `python setup.py install` in the main directory.

## Support and contributing

For bug reports and enhancement requests use the [Github issue tool](http://github.com/qimmuno/metaclonotypist/issues/new), or (even better!) open a [pull request](http://github.com/qimunno/metaclonotypist/pulls) with relevant changes. If you have any questions don't hesitate to contact us by email (qimmuno@gmail.com) or Twitter ([@andimscience](http://twitter.com/andimscience)).

You can run the testsuite by running `pytest` in the top-level directory.
"""

NAME = "metaclonotypist"
MAINTAINER = "Andreas Tiffeau-Mayer"
MAINTAINER_EMAIL = "andimscience@gmail.com"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://pyrepseq.readthedocs.io/"
DOWNLOAD_URL = "http://github.com/andim/pyrepseq"
LICENSE = "MIT"
AUTHOR = "Andreas Tiffeau-Mayer"
AUTHOR_EMAIL = "andimscience@gmail.com"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = find_packages()
PACKAGE_DATA = {"": [""]}
REQUIRES = [
    "numpy",
    "scipy",
    "pandas",
    "matplotlib",
    "pyrepseq",
]
