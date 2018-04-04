from setuptools import setup, find_packages
from codecs import open # To use a consistent encoding
from os.path import abspath, dirname, join, isfile, split
import os
import sys
import stat

here = abspath(dirname(__file__))

def get_shared_objects():
    r"""
    """
    cwd = dirname(abspath(__file__))
    sandy_dir = join(cwd, "sandy")
    list_so = [ fil for fil in os.listdir(sandy_dir) if '.so' in fil ]
    return list_so

def err():
    r"""
    Error function that exits with predifined message.
    """
    msg="THE INSTALLATION OF SANDY WAS NOT SUCCESSFUL"
    sys.exit(msg)


def makepdf():
    from subprocess import Popen, PIPE
    command = "make"
    args = [command, 'latexpdf']
    process = Popen(args, shell=False, cwd='docs', stdin=PIPE, stdout=PIPE,
                    stderr=PIPE)
    stdoutdata, stderrdata = process.communicate()
    if process.returncode != 0:
        sys.stdout.write(stdoutdata.decode())
        sys.stderr.write(stderrdata.decode())
        sys.stderr.write("ERROR IN THE COMPILATION OF SANDY'S DEVELOPMENT MANUAL")
        err()

#makepdf()


# Get the long description from the relevant file
#with open(join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

endf6_forlder = join("sandy", "endf6")
sampling_forlder = join("sandy", "sampling")
tests_folder = join("sandy", "data_tests")
#--------------------------------------------------------------------------
setup(
    name = 'sandy',
    version = '0.1',
    description = 'SANDY: sampling of nuclear data and uncertainty',
    #long_description = long_description,
    # The project's main homepage.
    # url='https://github.com/pypa/sampleproject',
    author = 'Luca Fiorito',
    author_email = 'lucafiorito.11@gmail.com',
#    license = '',
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        # 3 - Alpha
        # 4 - Beta
        # 5 - Production/Stable
        'Development Status :: 4 - Beta',
        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        # Pick your license as you wish (should match "license" above)
#        'License :: OSI Approved :: SCK-CEN',
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        ],
    keywords = 'uncertainty, nuclear data, covariance, sampling, sensitivity',
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
#    packages = ['sandy', njoy_folder],
    packages = find_packages(),

    # List run-time dependencies here. These will be installed by pip when your
    # project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/technical.html#install-requires-vs-requirements-files
    install_requires = ['numpy',
                        'scipy',
                        'pytest',
                        'matplotlib',
                        'pandas'],
    # If there are data files included in your packages that need to be
    # installed, specify them here.
    #package_data={ 'sandy' : ['*.so'], 'sandy.NJOY' : ['njoy.exe'],},
    include_package_data = True,

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
    'console_scripts': [
        'sandy=sandy.sampling.sampling:run'
        ],
    },
)
