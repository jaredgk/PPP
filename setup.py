from setuptools import setup
from Cython.Build import cythonize

requirements = [
      'pysam',
      'pandas',
      'pybedtools',
      'sklearn',
      'matplotlib'


]

setup(name="py-popgen",
      version="0.0.2",
      description="first setup file",
      include_package_data=True,
      packages=['pgpipe'],
      install_requires=requirements,
      scripts=['pgpipe/vcf_phase.py'],
      ext_modules = cythonize('pgpipe/test_cython.pyx'))

#packages=setuptools.find_packages() to automate package finding?
