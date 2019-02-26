import os
import sys
import gzip
import logging
from Bio import SeqIO

def check_fasta_gz (fasta_filename):

	try:
		fasta_file = gzip.open(fasta_filename)
		fasta_line = fasta_file.readline()
		fasta_file.close()

		if fasta_line[:1] == b'>':
			return True

		raise Exception('%s is not fasta formatted' % fasta_filename)

	except:

		return False

def check_fasta (fasta_filename):

	fasta_file = open(fasta_filename)
	fasta_line = fasta_file.readline()
	fasta_file.close()

	if fasta_line[:1] == b'>':
		return True

	raise Exception('%s is not fasta formatted' % fasta_filename)

def check_format(fasta_filename):
    """Checks header of given file for compression type


    Given a filename, opens file and reads first line to check if
    file has BGZF or GZIP header. May be extended to check for BCF format

    Parameters
    ----------
    fasta_filename : str
        Name of file to be checked

    Returns
    -------
    extension : str {'gzip', 'fasta', 'other'}
        File extension as indicated by header

    """

    # Check if the file is bgzipped or gunzipped, both of which should function
    if check_fasta_gz(fasta_filename):
   		return 'gzip'

    # Check if the file is in the fasta format
    if check_fasta(fasta_filename):
    	return 'fasta'
