"""
MDRestraintsGenerator
A framework for generating restraints for MD simulations
"""

# Add imports here
# from .restraints import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
