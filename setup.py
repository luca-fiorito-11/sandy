from setuptools import find_packages
from numpy.distutils.core import setup, Extension
import os


extensions = [
        Extension(name='rwf',
                  sources=[os.path.join(*['fortran', 'rwfortran.f'])]
                  ),
        ]
keywords = ['uncertainty', 'nuclear data', 'covariance', 'sampling', 'ENDF-6']
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
    keywords = ", ".join(keywords),
#    packages = find_packages(),
    data_files = [(x[0], list(map(lambda y: x[0]+'/'+y, x[2]))) for x in os.walk('sandy')],
    install_requires = open(requirements).read().splitlines(),
    zip_safe = False,
#    setup_requires=["pytest-runner",],
    tests_require=["pytest",],
    include_package_data = True,
    ext_modules = extensions,
    entry_points={
    'console_scripts': [
        'sandy=sandy.sampling:run',
        ],
    },
)
