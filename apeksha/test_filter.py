import unittest
import sys
import os
import logging
import filecmp
import shutil

import vcf_filter
import vcftools
import ruffus

class filterTest(unittest.TestCase):

    runpath= "/home/staff/asukeshkall/Downloads/Test/Run"
    if not os.path.exists(runpath):
       os.makedirs(runpath)
    

    def test_filter_run_1(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_1"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.bcf'))

        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.bcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

        

    def test_filter_run_2(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_2"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                    '--out', 'res'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.filter.log'))
        self.assertTrue(os.path.isfile('res.recode.bcf'))

        # Move the files to the new test folder
        shutil.move('res.filter.log', newpath)
        shutil.move('res.recode.bcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'res.filter.log')
        #self.addCleanup(os.remove, 'res.recode.bcf')

    def test_filter_run_3(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_3"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                    '--out-format', 'vcf'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.vcf'))

        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.vcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.vcf')


    def test_filter_run_4(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_4"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                    '--filter-max-alleles', '4'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.bcf'))
         
        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.bcf', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_filter_run_5(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_5"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                    '--filter-min-alleles', '1'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.bcf'))

        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.bcf', newpath)

         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_filter_run_6(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_6"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out', 'res',
                   '--out-format', 'vcf'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.filter.log'))
        self.assertTrue(os.path.isfile('res.recode.vcf'))
         
        # Move the files to the new test folder
        shutil.move('res.filter.log', newpath)
        shutil.move('res.recode.vcf', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'res.filter.log')
        #self.addCleanup(os.remove, 'res.recode.vcf')


    def test_filter_run_7(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_7"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out', 'res',
                   '--filter-max-alleles', '4'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.filter.log'))
        self.assertTrue(os.path.isfile('res.recode.bcf'))

        # Move the files to the new test folder
        shutil.move('res.filter.log', newpath)
        shutil.move('res.recode.bcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'res.filter.log')
        #self.addCleanup(os.remove, 'res.recode.bcf')

    def test_filter_run_8(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_8"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out', 'res',
                   '--filter-min-alleles', '1'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.filter.log'))
        self.assertTrue(os.path.isfile('res.recode.bcf'))

        # Move the files to the new test folder
        shutil.move('res.filter.log', newpath)
        shutil.move('res.recode.bcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'res.filter.log')
        #self.addCleanup(os.remove, 'res.recode.bcf')

    def test_filter_run_9(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_9"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out', 'res',
                   '--out-format','vcf',
                   '--filter-max-alleles', '4'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.filter.log'))
        self.assertTrue(os.path.isfile('res.recode.vcf'))

        # Move the files to the new test folder
        shutil.move('res.filter.log', newpath)
        shutil.move('res.recode.vcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'res.filter.log')
        #self.addCleanup(os.remove, 'res.recode.vcf')

    def test_filter_run_10(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_10"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out', 'res',
                   '--out-format','vcf',
                   '--filter-min-alleles', '1'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.filter.log'))
        self.assertTrue(os.path.isfile('res.recode.vcf'))

        # Move the files to the new test folder
        shutil.move('res.filter.log', newpath)
        shutil.move('res.recode.vcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'res.filter.log')
        #self.addCleanup(os.remove, 'res.recode.vcf')

    def test_filter_run_11(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_11"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out', 'res',
                   '--out-format','vcf',
                   '--filter-max-alleles', '4',
                   '--filter-min-alleles', '1'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.filter.log'))
        self.assertTrue(os.path.isfile('res.recode.vcf'))

        # Move the files to the new test folder
        shutil.move('res.filter.log', newpath)
        shutil.move('res.recode.vcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'res.filter.log')
        #self.addCleanup(os.remove, 'res.recode.vcf')

    def test_filter_run_12(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_12"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out-format','vcf',
                   '--filter-max-alleles', '4'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.vcf'))
         
        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.vcf', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.vcf')

    def test_filter_run_13(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_13"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out-format','vcf',
                   '--filter-max-alleles', '1'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.vcf'))

        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.vcf', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.vcf')

    def test_filter_run_14(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_14"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--out-format','vcf',
                   '--filter-max-alleles', '4',
                   '--filter-max-alleles', '1'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.vcf'))
         
        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.vcf', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.vcf')

    def test_filter_run_15(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/Run/test_15"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_filter.run(['merged_chr1_1000.vcf.gz',
                   '--filter-max-alleles', '4',
                   '--filter-max-alleles', '1'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.filter.log'))
        self.assertTrue(os.path.isfile('out.recode.bcf'))
         
        # Move the files to the new test folder
        shutil.move('out.filter.log', newpath)
        shutil.move('out.recode.bcf', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')



if __name__ == "__main__":
    unittest.main()
