import pysam
import argparse
import sys
from pgpipe.logging_module import initLogger
from pgpipe.model import read_model_file
from pgpipe.vcf_reader_func import VcfReader


def createParser():
    parser = argparse.ArgumentParser(description=("Remove individuals from a single VCF file"))
    parser.add_argument("--vcf",dest="vcfname",help=("VCF input filename"))
    parser.add_argument("--model-file",dest="modelname",help=("Model filename"))
    parser.add_argument("--model",dest="poptag",help=("Name of model if multiple in model file"))
    parser.add_argument("--out",dest="outname",default="-",help=("Name of VCF output file"))
    return parser






def filter_single_vcf(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs)
    popmodel = None
    if args.modelname is not None:
        popmodels = read_model_file(args.modelname)
        if len(popmodels) != 1:
            popmodel = popmodels[args.poptag]
        else:
            pp = list(popmodels.keys())
            popmodel = popmodels[pp[0]]
    vcf_in = VcfReader(args.vcfname,
                       popmodel=popmodel)
    sample_names = [l for l in vcf_in.reader.header.samples]
    bad_names = [False for l in sample_names]
    for rec in vcf_in.reader:
        #for i,samp in enumerate(rec.samples):
        for i in range(len(rec.samples)):
            samp = rec.samples[i]
            if samp.alleles[0] in [None,'N']:
                bad_names[i] = True
    final_names = []
    for i in range(len(sample_names)):
        if not bad_names[i]:
            final_names.append(sample_names[i])

    vcf_in.close()
    sys.stderr.write(str(len(final_names))+'\n')
    vcf_in = VcfReader(args.vcfname)
    vcf_in.reader.subset_samples(final_names)
    vcf_out = pysam.VariantFile(args.outname,'w',header=vcf_in.reader.header)
    rec_list = vcf_in.getRecordList()
    for rec in rec_list:
        vcf_out.write(rec)
    vcf_in.close()
    vcf_out.close()
    


if __name__ == "__main__":
    initLogger()
    filter_single_vcf(sys.argv[1:])