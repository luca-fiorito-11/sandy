from setuptools import setup, find_packages
from codecs import open # To use a consistent encoding
from os.path import abspath, dirname, join, isfile, split
import os
import sys
import stat


#def makepdf():
#    from subprocess import Popen, PIPE
#    command = "make"
#    args = [command, 'latexpdf']
#    process = Popen(args, shell=False, cwd='docs', stdin=PIPE, stdout=PIPE,
#                    stderr=PIPE)
#    stdoutdata, stderrdata = process.communicate()
#    if process.returncode != 0:
#        sys.stdout.write(stdoutdata.decode())
#        sys.stderr.write(stderrdata.decode())
#        sys.stderr.write("ERROR IN THE COMPILATION OF SANDY'S DEVELOPMENT MANUAL")
#        err()

#makepdf()


# Get the long description from the relevant file
#with open(join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

#--------------------------------------------------------------------------
setup(
    name = 'sandy',
    version = '0.1',
    description = 'SANDY: sampling of nuclear data and uncertainty',
    #long_description = long_description,
    # url='https://github.com/pypa/sampleproject',
    author = 'Luca Fiorito',
    author_email = 'lucafiorito.11@gmail.com',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        ],
    keywords = 'uncertainty, nuclear data, covariance, sampling, sensitivity',
    packages = find_packages(),
    install_requires = ['numpy',
                        'scipy',
                        'matplotlib',
                        'pytest',
                        'fortranformat>=0.2.5',
                        'pandas>=0.20'],
    include_package_data = True,
    entry_points={
    'console_scripts': [
        'sandy=sandy.sampling.sampling:run',
        'sandy_tests=sandy.sampling.tests:run',
        ],
    },
)
