import sys
import subprocess
import argparse


BGZIP_PATH='bgzip'
TABIX_PATH='tabix'

def bgzipFile(filename,force=False,renameout=None):
    args = [BGZIP_PATH]
    if force:
        args.append('-f')
    args.append(filename)
    if renameout is not None:
        args.append(['-c','>',renameout])
    print args
    bgzip_call = subprocess.Popen(args, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    bgzip_out, bgzip_err = bgzip_call.communicate()

def getFiletype(filename):
    try:
        file_ext = filename[-6:]
    except:
        sys.stderr.write("Filename %s is too short for valid extension\n"
                         % (filename))
    for ext in ['vcf','bed','sam','gff']:
        suff = ext+'.gz'
        if file_ext == suff:
            return ext
    raise Exception("Filename %s doesn't have a valid extension for tabix\n"
                    % (filename))


def tabixFile(filename,filetype='auto',force=True):
    args = [TABIX_PATH,filename]
    if filetype == 'auto':
        filetype = getFiletype(filename)
    if force:
        args.append('-f')
    print args
    tabix_call = subprocess.Popen(args, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    tabix_out, tabix_err = tabix_call.communicate()

def prepVcf(filename,force_b=False,force_t=True):
    bgzipFile(filename,force=force_b)
    compfilename = filename+'.gz'
    tabixFile(compfilename,force=force_t)


if __name__ == '__main__':
    prepVcf('example/mess.vcf')
    #tabixFile('example/mess.vcf.gz')
