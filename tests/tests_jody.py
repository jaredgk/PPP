import unittest
import filecmp
import sys
import os
import logging
import shutil
import zipfile
import gzip
import tempfile
import hashlib

# Import scripts to test

import pgpipe.vcf_bed_to_seq as vcf_bed_to_seq
import pgpipe.vcf_to_sfs as vcf_to_sfs
import pgpipe.vcf_to_dadi as vcf_to_dadi
import pgpipe.vcf_to_treemix as vcf_to_treemix
import pgpipe.vcf_to_gphocs as vcf_to_gphocs
import pgpipe.vcf_to_fastsimcoal as vcf_to_fastsimcoal


# Used to compare two files, returns bool
def file_comp(test_output, expected_output):
    return filecmp.cmp(test_output, expected_output)

# Used to compare two gz files, uses file_comp and a temp dir
def gz_file_comp (test_output, expected_output, tmp_dir):

    # Create the tmp paths
    tmp_test_path = os.path.join(tmp_dir, 'Test')
    tmp_expected_path = os.path.join(tmp_dir, 'Expected')

    # Create test tmp directories, if needed
    if not os.path.exists(tmp_test_path):
        os.makedirs(tmp_test_path)

    # Create test expected directories, if needed
    if not os.path.exists(tmp_expected_path):
        os.makedirs(tmp_expected_path)

    # Assign the tmp output files
    tmp_test_output = os.path.join(tmp_test_path, os.path.basename(test_output))
    tmp_expected_output = os.path.join(tmp_expected_path, os.path.basename(expected_output))

    # Open the gzip file
    with gzip.open(test_output, 'rb') as test_file:

        # Open the gunzip file
        with open(tmp_test_output, 'wb') as tmp_test_file:
            
            # Copy the file
            shutil.copyfileobj(test_file, tmp_test_file)

    # Open the gzip file
    with gzip.open(expected_output, 'rb') as expected_file:

        # Open the gunzip file
        with open(tmp_expected_output, 'wb') as tmp_expected_file:
            
            # Copy the file
            shutil.copyfileobj(expected_file, tmp_expected_file)

    # Check if the files have the same content
    file_compare_results = file_comp(tmp_test_output, tmp_expected_output)

    # Remove the tmp dirs
    shutil.rmtree(tmp_test_path)
    shutil.rmtree(tmp_expected_path)

    # Return the results
    return file_compare_results

class jh_function_tests (unittest.TestCase):

    # Assign the script directory
    @classmethod
    def setUpClass(cls):
        # Create a temporary directory
        cls.test_dir = tempfile.mkdtemp()
        cls.script_dir = os.path.dirname(os.path.realpath(__file__))
        cls.indir = os.path.join(cls.script_dir,'input')
        cls.expdir = os.path.join(cls.script_dir,'exp_output')
        
    @classmethod
    def tearDownClass (cls):

        # Remove the test directory after the tests
        shutil.rmtree(cls.test_dir)
        
    

    def test_vcf_bed_to_seq1 (self):
        vcf_bed_to_seq.run(vcf=os.path.join(self.indir,'pan_example.vcf.gz'),
                            fasta_reference=os.path.join(self.indir,"pan_example_ref.fa"),
                            model_file=os.path.join(self.indir,"panmodels.model"),
                            modelname = "4Pop",
                            region = "21:4431001-4499000",
                            out=os.path.join(self.expdir,'ppp_test_temp.out'))

       
        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_bed_to_seq_test.out')))

        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))


    def test_vcf_to_sfs1 (self):
        vcf_to_sfs.run(vcf=os.path.join(self.indir,"pan_example2.vcf.gz"),
               model_file=os.path.join(self.indir,"panmodels.model"),
               modelname = '4Pop',
               downsamplesizes = ['3','3','3','4'],
               folded = True,
               outgroup_fasta=os.path.join(self.indir,"chr22_pan_example2_ref.fa"),
               out=os.path.join(self.expdir,"ppp_test_temp.out"))

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_to_sfs_test1.txt')))

        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))

    def test_vcf_to_sfs2 (self):
        vcf_to_sfs.run(vcf=os.path.join(self.indir,"pan_example.vcf.gz"),
               model_file=os.path.join(self.indir,"panmodels.model"),modelname='5Pop',
               downsamplesizes=['3','3','3','4','2'],folded=True,
               outgroup_fasta=os.path.join(self.indir,"pan_example_ref.fa"),
               out=os.path.join(self.expdir,"ppp_test_temp.out"))
        
        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_to_sfs_test2.txt')))

        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))


    def test_vcf_to_gphocs1 (self):
        vcf_to_gphocs.run(vcf=os.path.join(self.indir,'pan_example.vcf.gz'),
                          reference=os.path.join(self.indir,'pan_example_ref.fa'),
                          model_file=os.path.join(self.indir,'panmodels.model'),
                          modelname = "4Pop",
                          bed_file=os.path.join(self.indir,'pan_example_regions.bed'),
                          out=os.path.join(self.expdir,'ppp_test_temp.out'))
        
        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_gphocs_test.out')))

        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))
    def test_vcf_to_gphocs2 (self):


        vcf_to_gphocs.run(vcf=os.path.join(self.indir,'pan_example2.vcf.gz'),
                          reference=os.path.join(self.indir,"chr22_pan_example2_ref.fa"),
                          model_file=os.path.join(self.indir,"panmodels.model"),
                          modelname="4Pop",
                          bed_file=os.path.join(self.indir,"pan_test_regions.bed"),
                          out=os.path.join(self.expdir,'ppp_test_temp.out'))
##        ,
##                          diploid=True)
        
        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_gphocs_test2.out')))

        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))      

    def test_vcf_to_treemix1(self):
        vcf_to_treemix.run(vcf = os.path.join(self.indir,'pan_example.vcf.gz'),
                        model_file = os.path.join(self.indir,'panmodels.model'),
                        modelname = '4Pop',
                        out = os.path.join(self.expdir,'ppp_test_temp.out'),
                        bed_file = os.path.join(self.indir,'pan_example_regions.bed'),
                        kblock = 1000)


        # Confirm that the output is what is expected
        self.assertTrue(gz_file_comp(os.path.join(self.expdir,'ppp_test_temp.out.gz'),os.path.join(self.expdir, 'vcf_treemixtest1.gz'), self.test_dir))

        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out.gz'))
        
    def test_vcf_to_treemix2(self):

        vcf_to_treemix.run(vcf = os.path.join(self.indir,'pan_example.vcf.gz'),
                        model_file = os.path.join(self.indir,'panmodels.model'),
                        modelname = '4Pop',
                        out = os.path.join(self.expdir,'ppp_test_temp.out'))
                           

        # Confirm that the output is what is expected
        self.assertTrue(gz_file_comp(os.path.join(self.expdir,'ppp_test_temp.out.gz'), os.path.join(self.expdir,'vcf_treemixtest2.gz'), self.test_dir))

        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out.gz'))        


    def test_vcf_to_fastsimcoal1(self):

        vcf_to_fastsimcoal.run(vcf = os.path.join(self.indir,"pan_example.vcf.gz"),
                    model_file = os.path.join(self.indir,"panmodels.model"),
                    modelname = '3Pop',
                    dim = ['1','2','m'],
                    basename =os.path.join(self.expdir,'ppp_test_temp.out'))

        # must extract archives to check that files match
        pz = zipfile.ZipFile(os.path.join(self.expdir,'ppp_test_temp.out.zip'),mode='r')
        vz = zipfile.ZipFile(os.path.join(self.expdir,'vcf_fsc1.zip'),mode='r')
        pznl = pz.namelist()
        vznl = vz.namelist()
        assert len(pznl) == len(vznl)
        for fi in range(len(pznl)):
            pz.extract(pznl[fi])
            vz.extract(vznl[fi])
        # Confirm that the output is what is expected
            self.assertTrue(file_comp(pznl[fi],vznl[fi]))
        # remove extracted files 
            self.addCleanup(os.remove, pznl[fi])
            self.addCleanup(os.remove, vznl[fi])          
      

    def test_vcf_to_fastsimcoal2(self):
        vcf_to_fastsimcoal.run(vcf = os.path.join(self.indir,"pan_example2.vcf.gz"),
                model_file = os.path.join(self.indir,"panmodels.model"),
                modelname = '5Pop',
                downsamplesizes = ['3','3','3','4','2'],
                basename = os.path.join(self.expdir,'ppp_test_temp.out'),
                folded = True,
                dim = ['1','2','m'],
                outgroup_fasta = os.path.join(self.indir,"chr22_pan_example2_ref.fa"))

        # must extract archives to check that files match 
        pz = zipfile.ZipFile(os.path.join(self.expdir,'ppp_test_temp.out.zip'),mode='r')
        vz = zipfile.ZipFile(os.path.join(self.expdir,'vcf_fsc2.zip'),mode='r')
        pznl = pz.namelist()
        vznl = vz.namelist()
        assert len(pznl) == len(vznl)
        for fi in range(len(pznl)):
            pz.extract(pznl[fi])
            vz.extract(vznl[fi])
        # Confirm that the output is what is expected
            self.assertTrue(file_comp(pznl[fi],vznl[fi]))
        # remove extracted files 
            self.addCleanup(os.remove, pznl[fi])
            self.addCleanup(os.remove, vznl[fi])  


    def test_vcf_to_dadi1(self):
        vcf_to_dadi.run(vcf = os.path.join(self.indir,'pan_example.vcf.gz'),
                model_file = os.path.join(self.indir,"panmodels.model"),
                modelname = "4Pop",
                out = os.path.join(self.expdir,'ppp_test_temp.out'),
                comment = 'testing bedfile',
                bed_file = os.path.join(self.indir,"pan_example_regions.bed"))
        
        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_dadisnp_bedfile_test.out')))
        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))

                        
    def test_vcf_to_dadi2(self):
        vcf_to_dadi.run(vcf = os.path.join(self.indir,"pan_example2.vcf.gz"),
                model_file = os.path.join(self.indir,"panmodels.model"),
                modelname = "4Pop",
                out = os.path.join(self.expdir,'ppp_test_temp.out'),
                comment = 'testing comment')
                        
        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_dadisnp_test.out')))
        # Remove the ouput
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))
                        
    def test_vcf_to_dadi3(self):
        vcf_to_dadi.run(vcf = os.path.join(self.indir,"pan_example2.vcf.gz"),
                        model_file = os.path.join(self.indir,"panmodels.model"),
                        modelname = "4Pop",
                        out  = os.path.join(self.expdir,'ppp_test_temp.out'),
                        comment = 'testing outgroup-fasta',
                        outgroup_fasta = os.path.join(self.indir,"chr22_pan_example2_ref.fa"))

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(os.path.join(self.expdir,'ppp_test_temp.out'),
                                  os.path.join(self.expdir,'vcf_dadisnp_fasta_test.out')))
        # Remove the ouput 
        self.addCleanup(os.remove, os.path.join(self.expdir,'ppp_test_temp.out'))

if __name__ == "__main__":
    unittest.main(verbosity = 2)



