import sys
import argparse
import gzip
import io
from vcf_reader_func import VcfReader
from model import Model, read_model_file

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--size',dest="window_size",type=int,default=1000)
    parser.add_argument('--zero-ho',dest="zero_ho",action="store_true",
                        help=("Output list in zero-based half-open format"))
    parser.add_argument('--zero-closed',dest="zero_closed",action="store_true",
                        help=("Output list in zero-based closed format"))
    parser.add_argument('--extend-regions',dest="extend",action="store_true",
                        help=("Regions will be cut at halfway point between valid and invalid site, instead of at site"))
    parser.add_argument('--missing-count',dest="missing_count",type=int,
                        help=("Count variants with more than (n) missing "
                        "samples as a missing site"),default="0")
    parser.add_argument('--vcf',dest="vcfname",type=str,
                        help=("Name of input VCF file, default is stdin"))
    parser.add_argument('--out',dest='outname',type=str,
                        help="Output filename (default is stdout)")
    parser.add_argument('--addchr',dest="addchr",action="store_true",
                        help=("Add 'chr' to start of chromosome name"))
    parser.add_argument('--removechr',dest="removechr",action="store_true",
                        help=("Remove 'chr' from start of chromosome name"))
    parser.add_argument('--pop',dest="popname",help="Model file if only using subset of inds in a file")
    parser.add_argument("--tbi", dest="tabix_index", help="Path to bgzipped file's index if name doesn't match VCF file")
    return parser

def getl(f,compressed=False):
    l = f.readline()
    if compressed:
        return l.decode()
    return l

def fixChromName(chrom,addchr,removechr):
    if addchr and chrom[3:] != 'chr':
        return 'chr'+chrom
    if removechr and chrom[3:] == 'chr':
        return chrom[3:]
    return chrom

def outputLine(chrom,start_pos,end_pos,args):
    sp = start_pos
    ep = end_pos
    if args.zero_ho:
        sp -= 1
    elif args.zero_closed:
        sp -= 1
        ep -= 1
    chrom = fixChromName(chrom,args.addchr,args.removechr)
    return chrom+'\t'+str(sp)+'\t'+str(ep)+'\n'

def regionsWithData(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs)
    compressed_input = False
    #if args.vcfname is None:
    #    instream = sys.stdin
    #elif args.vcfname[-3:] == '.gz':
    #    compressed_input = True
    #    instream = gzip.open(args.vcfname,'r')
   # else:
   #     instream = open(args.vcfname,'r')


    if args.vcfname is None:
        vcfname = '-'
    else:
        vcfname = args.vcfname

    popmodel = None
    if args.popname is not None:
        popmodels = read_model_file(args.modelname)
        pp = list(popmodels.keys())
        popmodel = popmodels[pp[0]]

    vcf_in = VcfReader(vcfname,
                       popmodel=popmodel,
                       index=args.tabix_index)

    if args.outname is not None:
        outstream = open(args.outname,'w')
    else:
        outstream = sys.stdout

    #sample_count = 20
    #missing_data_count = [0 for i in range(sample_count+1)]
    sample_count = 0
    indel_count = 0

    #missing_data_per_indiv = [0 for i in range(sample_count)]

    prev_miss_pos = 0
    prev_full_pos = 0
    prev_pos = 0
    prev_chrom = ''

    in_section = False
    #line = getl(instream,compressed_input)
    #while len(line.strip()) > 0:
    #for line in instream:
    for record in vcf_in.reader:
        #if line[0] == '#':
        #    line = getl(instream,compressed_input)
        #    continue
        #la = line.strip().split()
        missing_for_site = False
        missing_at_site = 0
        #for i in range(9,len(la)):
        for i in range(len(record.samples)):
            #geno = la[i]
            #if '.' in geno:
            if record.samples[i].alleles[0] in [None,'N']:
                missing_at_site += 1
                if missing_at_site > args.missing_count:
                    missing_for_site = True
                    break
        #cur_pos = int(la[1])
        cur_pos = record.pos
        #if prev_pos == 0 or prev_chrom != la[0]:
        if prev_pos == 0 or prev_chrom != record.chrom:
            if prev_pos != 0:
                start_pos = ((prev_miss_pos+prev_full_pos+1)//2 if args.extend else prev_full_pos)
                end_pos = prev_pos
                outstream.write(outputLine(prev_chrom,start_pos,end_pos,args))
            prev_miss_pos = cur_pos
            prev_full_pos = cur_pos
            #prev_chrom = la[0]
            prev_chrom = record.chrom
            if not missing_for_site:
                in_section = True
            prev_pos = cur_pos
            continue
        if missing_for_site:
            if in_section:
                in_section = False
                if args.extend:
                    start_pos = (prev_miss_pos+prev_full_pos+1)//2
                    end_pos = (prev_pos+cur_pos)//2
                else:
                    start_pos = prev_full_pos
                    end_pos = prev_pos
                diff = end_pos - start_pos
                if diff > args.window_size:
                    #if args.zero_ho:
                    #    start_pos -= 1
                    #elif args.zero_closed:
                    #    start_pos -= 1
                    #    end_pos -= 1
                    #chrom = fixChromName(la[0],args.addchr,args.removechr)
                    #outstream.write(chrom+'\t'+str(start_pos)+'\t'+str(end_pos)+'\n')
                    outstream.write(outputLine(record.chrom,start_pos,end_pos,args))
                    #sys.stdout.write(la[0]+'\t'+str(prev_full_pos)+'\t'+str(prev_pos)+'\n')
                    #print (la[0],prev_miss_pos,prev_full_pos,prev_pos,cur_pos,diff)
        else:
            if not in_section:
                prev_miss_pos = prev_pos
                prev_full_pos = cur_pos
                in_section = True
        prev_pos = cur_pos
        #line = getl(instream,compressed_input)
    if in_section:
        if args.extend:
            start_pos = (prev_miss_pos+prev_full_pos+1)//2
        else:
            start_pos = prev_full_pos
        end_pos = prev_pos
        outstream.write(outputLine(prev_chrom,start_pos,end_pos,args))
    try:
        outstream.close()
    except:
        pass

if __name__ == '__main__':
    regionsWithData(sys.argv[1:])
