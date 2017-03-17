import sys
from ruffus import *
import logging
from logging_module import formatLogger, switchLogger
from vcf_ref_to_seq import vcf_to_seq


start_list = open(str(sys.argv[1]),'r')
starting_files = [line.strip() for line in start_list]
#formatLogger()
#starting_files = ['example/chr11.subsamples.vcf.gz','example/chr11.vcf.gz']

#fmt_def = "%(asctime)s - %(levelname)s: %(message)s"
#fmtr = logging.Formatter(fmt=fmt_def)
#logging.basicConfig(filename='example/chr11.pipe.log',filemode="w",level=logging.INFO,format=fmt_def)

def exp_handler(type,val,tb):
    logging.exception("Error: %s" % (str(val)))

sys.excepthook = exp_handler

@transform(starting_files,
           suffix(".vcf.gz"),
           ".fasta")
def run_vcf_to_seq(vcf_in,seq_out):
    pref = vcf_in[0:-1*len(".vcf.gz")]
    switchLogger(pref)
    args = ['--vcf', vcf_in, '--ref', 'example/human_g1k_chr11.fasta',
            '--gr', 'example/snp_region.txt']
    vcf_to_seq(args)


@transform(run_vcf_to_seq,
           suffix(".fasta"),
           ".part")
def get_part(f_in,f_out):
    f = open(f_in,'r')
    o = open(f_out,'w')
    for line in f.readlines():
        o.write(line[0]+'\n')


def mainpipe(sysargs):
    pipeline_run()


if __name__ == '__main__':
    mainpipe(sys.argv[1:])
