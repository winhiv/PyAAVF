from setuptools import setup

from setuptools import find_packages, setup

DEPENDENCIES = ['setuptools']

setup(
    name='PyAAVF',
    version='0.1.0',
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    url='https://github.com/winhiv/PyAAVF.git',
    license='Apache License, Version 2.0',
    author='Matthew Fogel',
    author_email='matthew.fogel@canada.ca',
    install_requires=DEPENDENCIES,
    tests_require=(),
    description='An Amino Acid Variant Format parser for Python.',
    entry_points='',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ]
)
