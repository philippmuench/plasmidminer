mkdir CGAT
cd CGAT
wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
tar xvfz virtualenv-1.10.1.tar.gz
rm virtualenv-1.10.1.tar.gz
cd virtualenv-1.10.1
python virtualenv.py cgat-venv
source cgat-venv/bin/activate

# Install Python prerequisites
pip install cython

