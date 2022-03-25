from setuptools import setup
import os

defs = os.path.join('src')

with open('README.md', 'r') as f:
    long_description = f.read()

setup(name='XtraCrysPy',
      version='0.7',
      description='A Python tool for visualizing atomic systems and properties of condensed matter.',
      author='Frank T. Cerasoli',
      author_email='ftcerasoli@ucdavis.edu',
      platforms='Unix',
      url='https://github.com/sassafrass6/XtraCrysPy',
      packages=['XtraCrysPy'],
      package_dir={'XtraCrysPy':'src'},
      install_requires=['fury', 'matplotlib', 'numpy'],
      python_requires='>=3.8'
)
