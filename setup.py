from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
# https://python-packaging.readthedocs.io/en/latest/minimal.html
# For a fuller example see: https://github.com/CGATOxford/UMI-tools/blob/master/setup.py
# Or: https://github.com/CGATOxford/cgat/blob/master/setup.py

# See also this example: https://github.com/pypa/sampleproject/blob/master/setup.py

# TO DO: update with further options such as include README.rst and others when ready

# TO DO: to add tests see https://python-packaging.readthedocs.io/en/latest/testing.html

# TO DO: pass variable names for project and author when executing project_quickstart.py

from future import standard_library
standard_library.install_aliases()
from setuptools import setup

setup(name='BEST-D molecular analysis code',
      version='0.2',
      description='Code used in the BEST-D molecular analysis paper',
      url='https://github.com/AntonioJBT/BEST_D',
      author='Antonio J Berlanga-Taylor',
      author_email='a.berlanga at imperial.ac.uk',
      license='GPL-3.0',
#      packages=['funniest'],
#      install_requires=[
#            'cgat',
#            'CGATPipelines',
#      ],
      zip_safe=False
     )
