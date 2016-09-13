import os
from setuptools import find_packages
from setuptools import setup
import re


VERSIONFILE=os.path.join('plasmidminer', '_version.py')
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else: 
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

long_description = open('README.md', 'r').read()

setup(name='plasmidminer',
        version = verstr,
        description='plasmidminer',
        long_description = long_description,
        url = 'http://github.com/philippmuench/plasmidminer',
        author='Philipp C. Münch',
        author_email='philipp.muench@helmholtz-hzi.de',
        license='GNU General Public License, version 3 (GPL-3.0)',
        packages= ['plasmidminer'],
        include_package_data = True,
        scripts = ['plasmidminer/traitar', 'plasmidminer/merge_preds.py', 'plasmidminer/heatmap.py', 'plasmidminer/domtblout2gene_generic.py', 'traitar/predict.py', 'traitar/hmmer2filtered_best.py', 'traitar/hmm2gff.py'],
        zip_safe=False,
install_requires = ["matplotlib >= 1.3.1", "numpy >= 1.6", "scipy >= 0.13.3", "CGAT >=0.2.5"])

