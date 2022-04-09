from setuptools import setup
import os

defs = os.path.join('src')

long_description = open('./README.md', 'r').read()

setup(name='XtraCrysPy',
      version='0.7.4',
      description='A Python tool for visualizing atomic systems and properties of condensed matter.',
      long_description=long_description,
      author='Frank T. Cerasoli',
      author_email='ftcerasoli@ucdavis.edu',
      platforms='Unix',
      url='https://github.com/sassafrass6/XtraCrysPy',
      packages=['XtraCrysPy'],
      package_dir={'XtraCrysPy':'src'},
      install_requires=['fury', 'matplotlib', 'numpy'],
      python_requires='>=3.8'
)
