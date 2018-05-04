from setuptools import setup

DEPENDENCIES = ['setuptools']

setup(
    name='PyAAVF',
    version='0.1.0',
    packages=('PyAAVF',),
    url='https://github.com/winhiv/PyAAVF',
    license='MIT',
    author='Matthew Fogel',
    author_email='matthew.fogel@canada.ca',
    install_requires=DEPENDENCIES,
    tests_require=(),
    description='An Amino Acid Variant Format parser for Python.',
    entry_points=(
        'PyAAVF.filters': [
            'site_quality = PyAAVF.filters:SiteQuality',
            'vgq = PyAAVF.filters:AAVariantGenotypeQuality',
            'eb = PyAAVF.filters:ErrorBiasFilter',
            'dps = PyAAVF.filters:DepthPerSample',
            'avg-dps = PyAAVF.filters:AvgDepthPerSample',
            'snp-only = PyAAVF.filters.SnpOnly',
        ]
    )
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy')]
