"""
MDRestraintsGenerator
A framework for generating restraints for MD simulations
"""
from setuptools import setup, find_packages
import versioneer

with open('README.md', 'r') as handle:
    long_description = handle.read()

setup(
    # Self-descriptive entries which should always be present
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    name='MDRestraintsGenerator',
    author='Irfan Alibay',
    author_email='irfan.alibay@bioch.ox.ac.uk',
    description='Enabling the use of restraints in alchemical simulations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='LGPLv3',
    url='https://github.com/IAlibay/MDRestraintsGenerator/',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.8',
    install_requires=[
        'MDAnalysis>1.0.0',
        'numpy',
        'scipy<1.8',
        'matplotlib',
    ],
)
