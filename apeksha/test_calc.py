import unittest
import sys
import os
import logging
import filecmp
import shutil

import vcf_calc
import vcftools
import ruffus

class filterTest(unittest.TestCase):

    runpath= "/home/staff/asukeshkall/Downloads/Test/run_calc"
    if not os.path.exists(runpath):
       os.makedirs(runpath)
    
    def test_calc_run_1(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_1"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_2(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_2"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('res.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.weir.fst', newpath)
        shutil.move('res.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_3(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_3"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'weir-fst',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.weir.fst'))
        self.assertTrue(os.path.isfile('out.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.weir.fst', newpath)
        shutil.move('out.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_4(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_4"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_5(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_5"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'TajimaD',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.Tajima.D'))
        self.assertTrue(os.path.isfile('out.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('out.Tajima.D', newpath)
        shutil.move('out.Tajima.D.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_6(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_6"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'pi',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.pi'))
        self.assertTrue(os.path.isfile('out.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.pi', newpath)
        shutil.move('out.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_7(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_7"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'freq',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.frq'))
        self.assertTrue(os.path.isfile('out.frq.log'))

        # Move the files to the new test folder
        shutil.move('out.frq', newpath)
        shutil.move('out.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_8(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_8"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'het',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.het'))
        self.assertTrue(os.path.isfile('out.het.log'))

        # Move the files to the new test folder
        shutil.move('out.het', newpath)
        shutil.move('out.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_9(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_9"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_10(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_10"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_11(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_11"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'weir-fst',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.weir.fst'))
        self.assertTrue(os.path.isfile('res.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.weir.fst', newpath)
        shutil.move('res.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_12(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_12"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('res.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.weir.fst', newpath)
        shutil.move('res.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_13(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_13"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'TajimaD',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.Tajima.D'))
        self.assertTrue(os.path.isfile('res.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('res.Tajima.D', newpath)
        shutil.move('res.Tajima.D.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_14(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_14"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'pi',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.pi'))
        self.assertTrue(os.path.isfile('res.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.pi', newpath)
        shutil.move('res.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_15(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_15"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'freq',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.frq'))
        self.assertTrue(os.path.isfile('res.frq.log'))

        # Move the files to the new test folder
        shutil.move('res.frq', newpath)
        shutil.move('res.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_16(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_16"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'het',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.het'))
        self.assertTrue(os.path.isfile('res.het.log'))

        # Move the files to the new test folder
        shutil.move('res.het', newpath)
        shutil.move('res.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_17(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_17"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('res.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.weir.fst', newpath)
        shutil.move('res.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_18(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_18"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('res.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.weir.fst', newpath)
        shutil.move('res.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_19(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_19"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'weir-fst',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.weir.fst'))
        self.assertTrue(os.path.isfile('res.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.weir.fst', newpath)
        shutil.move('res.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_20(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_20"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('res.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.weir.fst', newpath)
        shutil.move('res.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_21(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_21"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'TajimaD',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.Tajima.D'))
        self.assertTrue(os.path.isfile('res.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('res.Tajima.D', newpath)
        shutil.move('res.Tajima.D.log', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_22(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_22"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'pi',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.pi'))
        self.assertTrue(os.path.isfile('res.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.pi', newpath)
        shutil.move('res.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_23(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_23"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'freq',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.frq'))
        self.assertTrue(os.path.isfile('res.frq.log'))

        # Move the files to the new test folder
        shutil.move('res.frq', newpath)
        shutil.move('res.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_24(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_24"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'het',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.het'))
        self.assertTrue(os.path.isfile('res.het.log'))

        # Move the files to the new test folder
        shutil.move('res.het', newpath)
        shutil.move('res.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_25(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_25"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'weir-fst',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.weir.fst'))
        self.assertTrue(os.path.isfile('res.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.weir.fst', newpath)
        shutil.move('res.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_26(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_26"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('res.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.weir.fst', newpath)
        shutil.move('res.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_27(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_27"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'TajimaD',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.Tajima.D'))
        self.assertTrue(os.path.isfile('res.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('res.Tajima.D', newpath)
        shutil.move('res.Tajima.D.log', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_28(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_28"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'pi',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.pi'))
        self.assertTrue(os.path.isfile('res.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.pi', newpath)
        shutil.move('res.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_29(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_29"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'freq',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.frq'))
        self.assertTrue(os.path.isfile('res.frq.log'))

        # Move the files to the new test folder
        shutil.move('res.frq', newpath)
        shutil.move('res.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_30(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_30"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'het',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.het'))
        self.assertTrue(os.path.isfile('res.het.log'))

        # Move the files to the new test folder
        shutil.move('res.het', newpath)
        shutil.move('res.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_31(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_31"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'weir-fst',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.weir.fst'))
        self.assertTrue(os.path.isfile('res.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.weir.fst', newpath)
        shutil.move('res.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_32(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_32"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('res.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.weir.fst', newpath)
        shutil.move('res.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_33(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_33"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'TajimaD',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.Tajima.D'))
        self.assertTrue(os.path.isfile('res.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('res.Tajima.D', newpath)
        shutil.move('res.Tajima.D.log', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_34(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_34"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'pi',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.windowed.pi'))
        self.assertTrue(os.path.isfile('res.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('res.windowed.pi', newpath)
        shutil.move('res.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_35(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_35"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'freq',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.frq'))
        self.assertTrue(os.path.isfile('res.frq.log'))

        # Move the files to the new test folder
        shutil.move('res.frq', newpath)
        shutil.move('res.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_36(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_36"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--out', 'res',
                      '--calc-statistic', 'het',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('res.het'))
        self.assertTrue(os.path.isfile('res.het.log'))

        # Move the files to the new test folder
        shutil.move('res.het', newpath)
        shutil.move('res.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')


    def test_calc_run_37(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_37"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'weir-fst',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.weir.fst'))
        self.assertTrue(os.path.isfile('out.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.weir.fst', newpath)
        shutil.move('out.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_38(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_38"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_39(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_39"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'TajimaD',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.Tajima.D'))
        self.assertTrue(os.path.isfile('out.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('out.Tajima.D', newpath)
        shutil.move('out.Tajima.D.log', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_40(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_40"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'pi',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.pi'))
        self.assertTrue(os.path.isfile('out.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.pi', newpath)
        shutil.move('out.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_41(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_41"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'freq',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.frq'))
        self.assertTrue(os.path.isfile('out.frq.log'))

        # Move the files to the new test folder
        shutil.move('out.frq', newpath)
        shutil.move('out.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_42(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_42"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'het',
                      '--statistic-window-size', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.het'))
        self.assertTrue(os.path.isfile('out.het.log'))

        # Move the files to the new test folder
        shutil.move('out.het', newpath)
        shutil.move('out.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_43(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_43"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'weir-fst',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.weir.fst'))
        self.assertTrue(os.path.isfile('out.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.weir.fst', newpath)
        shutil.move('out.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_44(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_44"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_45(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_45"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'TajimaD',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.Tajima.D'))
        self.assertTrue(os.path.isfile('out.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('out.Tajima.D', newpath)
        shutil.move('out.Tajima.D.log', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_46(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_46"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'pi',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.pi'))
        self.assertTrue(os.path.isfile('out.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.pi', newpath)
        shutil.move('out.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_47(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_47"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'freq',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.frq'))
        self.assertTrue(os.path.isfile('out.frq.log'))

        # Move the files to the new test folder
        shutil.move('out.frq', newpath)
        shutil.move('out.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_48(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_48"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'het',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.het'))
        self.assertTrue(os.path.isfile('out.het.log'))

        # Move the files to the new test folder
        shutil.move('out.het', newpath)
        shutil.move('out.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_49(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_49"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'weir-fst',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.weir.fst'))
        self.assertTrue(os.path.isfile('out.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.weir.fst', newpath)
        shutil.move('out.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_50(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_50"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_51(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_51"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'TajimaD',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.Tajima.D'))
        self.assertTrue(os.path.isfile('out.Tajima.D.log'))

        # Move the files to the new test folder
        shutil.move('out.Tajima.D', newpath)
        shutil.move('out.Tajima.D.log', newpath)

        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_52(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_52"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'pi',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.pi'))
        self.assertTrue(os.path.isfile('out.windowed.pi.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.pi', newpath)
        shutil.move('out.windowed.pi.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_53(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_53"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'freq',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.frq'))
        self.assertTrue(os.path.isfile('out.frq.log'))

        # Move the files to the new test folder
        shutil.move('out.frq', newpath)
        shutil.move('out.frq.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_54(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_54"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--calc-statistic', 'het',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.het'))
        self.assertTrue(os.path.isfile('out.het.log'))

        # Move the files to the new test folder
        shutil.move('out.het', newpath)
        shutil.move('out.het.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')

    def test_calc_run_55(self):
        newpath= "/home/staff/asukeshkall/Downloads/Test/run_calc/test_55"
        if not os.path.exists(newpath):
           os.makedirs(newpath)

        vcf_calc.run(['merged_chr1_1000.vcf.gz',
                      '--statistic-window-size', '2',
                      '--statistic-window-step', '2',
                      '--pop-file', 'Paniscus.txt',
                      '--pop-file', 'Troglodytes.txt'])

        # Confirm that the output is what is expected      
        self.assertTrue(os.path.isfile('out.windowed.weir.fst'))
        self.assertTrue(os.path.isfile('out.windowed.weir.fst.log'))

        # Move the files to the new test folder
        shutil.move('out.windowed.weir.fst', newpath)
        shutil.move('out.windowed.weir.fst.log', newpath)
         
        # Remove the ouput and log files created by the function
        #self.addCleanup(os.remove, 'out.filter.log')
        #self.addCleanup(os.remove, 'out.recode.bcf')





if __name__ == "__main__":
    unittest.main()
