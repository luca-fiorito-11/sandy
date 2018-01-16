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
    msg="THE INSTALLATION OF SANDY WAS NOT SUCCESFUL"
    sys.exit(msg)


def compile_fortran():
    r"""
    Run ``f2py`` from the shell to confirm that it is installed (it works 
    only in LINUX).
    
    If ``f2py`` is installed, then compile module ``sandy.rw_fortran.f`` into 
    a shared object (``.so``).
    """
    from numpy.f2py import compile
    from shutil import move
    f = open( join('sandy','rw_fortran.f'), 'r' )
    SOURCE = f.read()
    f.close()
    if sys.platform == 'win32':
        eargs = '--compiler=mingw32'
    else:
        eargs = ''
    compile(SOURCE.encode(), modulename = 'rw_fortran', extra_args = eargs, verbose = 0)
    for f in os.listdir():
        if 'rw_fortran' in f:
            if f in os.listdir(join(here,'sandy')):
                sys.stdout.write("FILE '{}' ALREADY EXISTS".format(f))
            else:
                move(f, join(here,'sandy'))
            return
    sys.stderr.write("COULD NOT COPY SRC='{}' TO DTS='{}'".format(njoy_source, njoy))
    err()


def get_njoy():
    """
    Get ``NJOY`` executable from keyboard input and copy it to correct folder 
    (``sandy.NJOY``) with new name ``njoy.exe``.
    Then make it executable.
    """
    from shutil import copy
    njoy_source = input("\nPROVIDE NJOY EXECUTABLE [e.g. path/to/njoy.exe] : ")
    if not isfile(njoy_source):
        sys.exit("NJOY EXECUTABLE FILE '{}' DOES NOT EXIST".format(njoy_source))
    njoy = join(here, join('sandy', join('NJOY', 'njoy.exe')))
    try:
        copy(njoy_source, njoy)
    except:
        sys.stderr.write("COULD NOT COPY SRC='{}' TO DTS='{}'".format(njoy_source, njoy))
        err()
    st = os.stat(njoy)
    os.chmod(njoy, st.st_mode | stat.S_IEXEC)

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
compile_fortran()
get_njoy()


# Get the long description from the relevant file
#with open(join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

tests_folder = join("sandy", "tests")
njoy_folder = join("sandy", "NJOY")
#--------------------------------------------------------------------------
setup(
    name = 'sandy',
    version = '0.1',
    description = 'SANDY: sampling of nuclear data and uncertainty',
    #long_description = long_description,
    # The project's main homepage.
    # url='https://github.com/pypa/sampleproject',
    author = 'Luca Fiorito',
    author_email = 'luca.fiorito@sckcen.be',
    license = 'SCK-CEN',
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
        'License :: OSI Approved :: SCK-CEN',
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
    install_requires = ['numpy', 'scipy', 'pytest', 'matplotlib'],

    # If there are data files included in your packages that need to be
    # installed, specify them here.
    #package_data={ 'sandy' : ['*.so'], 'sandy.NJOY' : ['njoy.exe'],},
    include_package_data = True,

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
    'console_scripts': [
        'sandy=sandy.run:preprocessing',
        'sandy_tests=sandy.test_modules:run_tests'
        ],
    },
)
