"""
    reads an ima3 input file
    generates a gzipped treemix file
    takes two arguments at the command line
        input file
        output file

    the output file is to be run using treemix's block option
    the data from a locus in the ima3 file is padded out to have a length of
    200 snps  (i.e. padded with uninformative snps)

    Then treemix is run with -k 200

"""
import sys
import gzip

def getVariantSites(seq_list):
    idx = []
    for i in range(len(seq_list[0])):
        first_base = seq_list[0][i]
        for j in range(1,len(seq_list)):
            cur_base = seq_list[j][i]
            if first_base != cur_base:
                idx.append(i)
                break
    return idx



def parseSeqs(seq_list, seq_count, padline,padcount):
    idx = getVariantSites(seq_list)
    var_list = []
    for i in idx:
        ref_count = [0 for i in range(len(seq_count))]
        alt_count = [0 for i in range(len(seq_count))]
        ref = seq_list[0][i]
        #ref_count[0] = 1
        cur_seq = 0
        for j in range(len(seq_count)):
            for k in range(seq_count[j]):
                if seq_list[cur_seq][i] == ref:
                    ref_count[j] += 1
                else:
                    alt_count[j] += 1
                cur_seq += 1
        var_line = ' '.join([str(a)+','+str(b) for a,b in zip(ref_count,alt_count)])
        var_list.append(var_line)
    if padcount < len(var_list):
        raise Exception("too many snps")
    while len(var_list) < padcount:
        var_list.append(padline)
    return var_list

def main(imfilename,treemixfilename):

    padcount = 200
    f = open(str(imfilename),'r')
    locus_list = []

    pop_list = []

    output_contents = []

    l = f.readline()
    l = f.readline()
    while l[0] == '#':
        l = f.readline()
    npops = l.split()[0]
    pop_list = f.readline().strip().split()
    print ("populations",pop_list)
    output_contents.append(' '.join(pop_list))
    temp = f.readline().strip()
    print ("number of loci:",temp)
    if temp.find('(') >= 0:
        treestring = temp
        num_loci = int(f.readline().strip())
    else:
        num_loci = int(temp)
    for i in range(num_loci):
        loci_data = f.readline().strip()
        loci_array = loci_data.split()
        seq_count = [int(ll) for ll in loci_array[1:1+len(pop_list)]]
        seq_total = sum(seq_count)
        padline = ""
        for c in seq_count:
            padline  = padline + str(c) + ",0 "
        seq_list = []
        for j in range(seq_total):
            seq_list.append(f.readline().strip().split()[1])
        sl = parseSeqs(seq_list, seq_count,padline[0:-1],padcount)
        for s in sl:
            output_contents.append(s)

    out_buffer = '\n'.join(output_contents)+'\n'
    if str(treemixfilename).endswith('gz'):
        with gzip.open(str(treemixfilename), 'wb') as out_f:
            out_f.write(out_buffer.encode())
    else:
        with gzip.open(str(treemixfilename)+".gz", 'wb') as out_f:
            out_f.write(out_buffer.encode())



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("reads an ima3 input file\ngenerates a gzipped treemix file\ntakes two arguments at the command line\n\tinput file\toutput file")
    main(sys.argv[1],sys.argv[2])

