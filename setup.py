"""
MDRestraintsGenerator
A framework for generating restraints for MD simulations
"""
from setuptools import setup
import versioneer


setup(
    # Self-descriptive entries which should always be present
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
