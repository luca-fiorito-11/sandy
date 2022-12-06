from setuptools import find_packages
from numpy.distutils.core import setup, Extension
import os
import sys
import logging

__author__ = "Luca Fiorito"

keywords = ['uncertainty', 'nuclear data', 'covariance', 'sampling', 'ENDF-6']
requirements = "requirements.txt"

setup(
      name='sandy',
      version='0.1',
      description='SANDY: sampling of nuclear data and uncertainty',
      url='https://github.com/luca-fiorito-11/sandy',
      author='Luca Fiorito',
      author_email='lucafiorito.11@gmail.com',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Programming Language :: Python :: 3',
          ],
      keywords=", ".join(keywords),
      # packages = find_packages(),
      data_files=[(x[0], list(map(lambda y: x[0]+'/'+y, x[2]))) for x in os.walk('sandy')],
      install_requires=open(requirements).read().splitlines(),
      zip_safe=False,
      # setup_requires=["pytest-runner",],
      tests_require=[
          "pytest",
          ],
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'sandy=sandy.sampling:run',
              ],
          },

      )

#import sandy
#
#try:
#    sandy.get_njoy()
#except sandy.SandyError:
#    logging.warning("env variable 'NJOY' is not assigned. SANDY might not behave as expected.")
