from numpy.distutils.core import setup
import os

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
      data_files=[(x[0], list(map(lambda y: x[0]+'/'+y, x[2]))) for x in os.walk('sandy')],
      install_requires=open(requirements).read().splitlines(),
      zip_safe=False,
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
