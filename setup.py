from distutils.core import setup

version = "0.1"

setup(name='plasmidminer',
        version = version,
        description='mining of plasmids from metagenomic samples',
        long_description = 'some text goes here',
        url = 'http://github.com/philippmuench/plasmidminer',
        author='Philipp C. Münch',
        author_email='philipp.muench@helmholtz-hzi.de',
        packages= ['plasmidminer'],
        include_package_data = True,
        scripts = ['plasmidminer/plasmidminer', 'plasmidminer/download.py', 'plasmidminer/features.py', 'plasmidminer/simulate_metagenomic_reads.py', 'traitar/simulate_reads.py'],
        zip_safe=False,
	install_requires = ["matplotlib >= 1.3.1", "numpy >= 1.6", "scipy >= 0.13.3", "CGAT >=0.2.5"]
)

