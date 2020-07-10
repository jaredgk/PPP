import io
import pgpipe
from setuptools import setup
#from Cython.Build import cythonize

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
               'pgpipe/vcf_four_gamete.py',
               'pgpipe/stat_sampler.py',
               'pgpipe/bed_utilities.py',
               'pgpipe/informative_loci_filter.py',
               'pgpipe/vcf_format_conversions.py',
               'pgpipe/vcf_to_ima.py',
               'pgpipe/admixture.py',
               'pgpipe/plink_linkage_disequilibrium.py',
               'pgpipe/ima3_wrapper.py',
               'pgpipe/eigenstrat_fstats.py',
               'pgpipe/vcf_utilities.py',
               'pgpipe/vcf_split_pysam.py',
               'pgpipe/vcf_to_fastsimcoal.py',
               'pgpipe/vcf_to_dadi.py',
               'pgpipe/vcf_to_treemix.py',
               'pgpipe/vcf_to_gphocs.py',
               'pgpipe/model_creator.py',
               'pgpipe/vcf_to_sfs.py',
               'pgpipe/vcf_bed_to_seq.py']

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
      python_requires=">=3.6.*, !=3.8.*")
#     , ext_modules = cythonize('pgpipe/test_cython.pyx'))

#packages=setuptools.find_packages() to automate package finding?
