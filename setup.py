import io
import pgpipe
from setuptools import setup
from Cython.Build import cythonize

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

requirements = ['Cython',
                'pysam',
                'pandas',
                'pybedtools',
                'scikit-learn',
                'matplotlib']

ppp_scripts = ['pgpipe/vcf_filter.py',
               'pgpipe/vcf_calc.py',
               'pgpipe/vcf_split.py',
               'pgpipe/vcf_phase.py',
               'pgpipe/four_gamete.py',
               'pgpipe/stat_sampler.py',
               'pgpipe/bed_utilities.py',
               'pgpipe/informative_loci_filter.py',
               'pgpipe/convert.py',
               'pgpipe/vcf_to_ima.py',
               'pgpipe/admixture.py',
               'pgpipe/plink_ld.py',
               'pgpipe/ima3_wrapper.py',
               'pgpipe/eigenstrat_fstats.py',
               'pgpipe/vcf_liftover.py',
               'pgpipe/vcf_utilities.py']

setup(name=pgpipe.__title__,
      version=pgpipe.__version__,
      project_urls={"Documentation": "https://ppp.readthedocs.io",
                    "Code": "https://github.com/jaredgk/PPP",
                    "Issue tracker": "https://github.com/jaredgk/PPP/issues"},
      license=pgpipe.__license__,
      url=pgpipe.__url__,
      author=pgpipe.__author__,
      author_email=pgpipe.__email__,
      maintainer="Jared Knoblauch, " \
                 "Andrew Webb",
      maintainer_email="tug41380@temple.edu",
      description=pgpipe.__summary__,
      long_description=readme,
      include_package_data=True,
      packages=['pgpipe'],
      install_requires=requirements,
      scripts=ppp_scripts,
      python_requires=">=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5.*",
      ext_modules = cythonize('pgpipe/test_cython.pyx'))

#packages=setuptools.find_packages() to automate package finding?
