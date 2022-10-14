from setuptools import setup
import os

defs = os.path.join('src')

long_description = open('./README.md', 'r').read()

setup(name='XtraCrysPy',
      version='1.0.2',
      description='A Python tool for visualizing atomic systems and properties of condensed matter.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Frank T. Cerasoli',
      author_email='ftcerasoli@ucdavis.edu',
      platforms='Unix',
      url='https://github.com/sassafrass6/XtraCrysPy',
      packages=['XtraCrysPy'],
      package_dir={'XtraCrysPy':'src'},
      install_requires=['ase', 'fury', 'numpy'],
      python_requires='>=3.8'
)
