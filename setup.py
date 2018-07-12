from setuptools import find_packages
from numpy.distutils.core import setup, Extension
import os


extensions = [
        Extension(name='rwf',
                  sources=[os.path.join(*['fortran', 'rwfortran.f'])]
                  ),
        ]
requirements = "requirements.txt"
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
    install_requires = open(requirements).read().splitlines(),
    zip_safe = False,
#    ['numpy',
#                        'scipy',
#                        'matplotlib',
#                        'pytest>=3.3',
#                        'fortranformat>=0.2.5',
#                        'bokeh>=0.12.10',
#                        'pandas>=0.20'],
#    setup_requires=["pytest-runner",],
    tests_require=["pytest",],
    include_package_data = True,
    ext_modules = extensions,
#    entry_points={
#    'console_scripts': [
#        'sandy=sandy.sampling.sampling:run',
#        'sandy_tests=sandy.sampling.tests:runtests',
#        'sandy_xs_plotter=sandy.sampling.plotter2:main',
#        'sandy_njoy=sandy.njoy.njoy:process_lib'
#        ],
#    },
)
