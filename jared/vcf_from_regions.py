import sys
import pysam
import argparse
import logging
from logging_module import individualFunctionLogger
from gene_region import Region, RegionList
import vcf_reader_func as vf



def createParser():
    parser.argparse.ArgumentParser(description=("Given a range or a file of "
                                "ranges and a VCF file, will generate one or "
                                "more VCF files with variants only from the "
                                "region(s) specified"))
    parser.add
