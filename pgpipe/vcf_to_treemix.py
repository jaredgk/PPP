import sys
import argparse
import pysam
import gzip

from pgpipe.vcf_reader_func import VcfReader
from pgpipe.model import Model, read_single_model

def createParser():
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument("--vcfs",nargs="+",help=("VCFs that represent loci for input"))
    parser.add_argument("--model-file",dest="model_file",help="Name of model file")
    parser.add_argument("--model",dest="model",help="Name of model in model file if more than one model is present")
    parser.add_argument("--out",dest="outname",help="Name of output treemix file")
    parser.add_argument("--pad",dest="pad",type=int,default=200,help=("Loci will be padded to this many sites (must be larger than locus with most variant sites"))
    parser.add_argument("--compress-out",dest="compress",action="store_true",help=("Will compress output, adding .gz to output name if not already specified"))

    return parser

def createVarString(rec,vcf_reader):
    variant_site = False
    outstr = ''
    for p in vcf_reader.popkeys.keys():
        refcount = 0
        altcount = 0
        for indiv_idx in vcf_reader.popkeys[p]:
            for hap in range(len(rec.samples[indiv_idx].alleles)):
                if rec.samples[indiv_idx].alleles[hap] != rec.ref:
                    altcount += 1
                    variant_site = True
                else:
                    refcount += 1
        outstr += str(refcount)+','+str(altcount)+' '
    return outstr.strip(),variant_site

def createInvarString(rec,vcf_reader):
    outstr = ''
    for p in vcf_reader.popkeys.keys():
        refcount = 0
        for indiv_idx in vcf_reader.popkeys[p]:
            refcount += len(rec.samples[indiv_idx].alleles)
        outstr += str(refcount)+',0 '
    return outstr.strip()

def createPopString(vcf_reader):
    return ' '.join(vcf_reader.popkeys.keys())

def writeToFile(outf,line,compress):
    if compress:
        line = line.encode()
    outf.write(line)

def vcf_treemix_convert(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs)
    popmodel = read_single_model(args.model_file,args.model)
    if args.compress:
        outfn = (args.out+'.gz' if args.out[-3:] != '.gz' else args.out)
        output_file = gzip.open(outfn,'wb')
    else:
        output_file = open(args.outname,'w')

    pop_out = False

    for vcfname in args.vcfs:
        vcf_reader = VcfReader(vcfname,popmodel=popmodel)
        rec_list = vcf_reader.getRecordList()
        vcf_reader.reader.close()
        var_count = 0
        if not pop_out:
            writeToFile(output_file,createPopString(vcf_reader)+'\n',args.compress)
            pop_out = True
        for rec in rec_list:
            sitestr,varsite = createVarString(rec,vcf_reader)
            if varsite:
                writeToFile(output_file,sitestr+'\n',args.compress)
                var_count += 1
        if var_count > args.pad:
            raise Exception("%s has %d loci, over max of %d"%(vcfname,var_count,args.pad))
        for i in range(var_count,args.pad):
            writeToFile(output_file,createInvarString(rec_list[0],vcf_reader)+'\n',args.compress)







if __name__ == "__main__":
    vcf_treemix_convert(sys.argv[1:])