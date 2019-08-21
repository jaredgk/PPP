import io
from pgpipe import __version__
from setuptools import setup
from Cython.Build import cythonize

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

author_str = "Andrew Webb, " \
             "Jared Knoblauch, " \
             "Nitesh Sabankar, " \
             "Apeksha Sukesh Kallur, " \
             "Jody Hey, " \
             "Arun Sethuraman"

email_str = "tug41380@temple.edu, " \
            "jaredknoblauch@gmail.com, " \
            "nsabankar@csusm.edu, " \
            "asukeshkall@csusm.edu, " \
            "hey@temple.edu, " \
            "asethuraman@csusm.edu"

requirements = ['pysam',
                'pandas',
                'pybedtools',
                'sklearn',
                'matplotlib']

ppp_scripts = ['pgpipe/vcf_filter.py',
               'pgpipe/vcf_calc.py',
               'pgpipe/vcf_split.py',
               'pgpipe/vcf_phase.py',
               'pgpipe/four_gamete.py',
               'pgpipe/stat_sampler.py',
               'pgpipe/bed_tasks.py',
               'pgpipe/informative_loci_filter.py',
               'pgpipe/convert.py',
               'pgpipe/vcf_to_ima.py',
               'pgpipe/admixture.py',
               'pgpipe/plink_ld.py',
               'pgpipe/ima3_wrapper.py',
               'pgpipe/eigenstrat_fstats.py',
               'pgpipe/vcf_liftover.py',
               'pgpipe/vcf_utilities.py']

setup(name="py-popgen",
      version=__version__,
      project_urls={"Documentation": "https://ppp.readthedocs.io",
                    "Code": "https://github.com/jaredgk/PPP",
                    "Issue tracker": "https://github.com/jaredgk/PPP/issues"},
      license="MIT",
      author=author_str,
      author_email=email_str,
      maintainer="Jared Knoblauch, " \
                 "Andrew Webb",
      maintainer_email="jaredknoblauch@gmail.com, " \
                       "tug41380@temple.edu",
      description="Software platform for facilitating population genomic analyses",
      long_description=readme,
      include_package_data=True,
      packages=['pgpipe'],
      install_requires=requirements,
      scripts=ppp_scripts,
      python_requires=">=2.7, >=3.6",
      ext_modules = cythonize('pgpipe/test_cython.pyx'))

#packages=setuptools.find_packages() to automate package finding?
