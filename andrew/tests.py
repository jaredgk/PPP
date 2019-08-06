import unittest
import filecmp
import sys
import os
import logging
import shutil
import numpy as np

# Import scripts to test
import pgpipe.vcftools as vcftools
import pgpipe.vcf_calc as vcf_calc
import pgpipe.vcf_sampler as vcf_sampler

# Used to compare two files, returns bool
def file_comp(test_output, expected_output):
    return filecmp.cmp(test_output, expected_output)

# Run tests for the functions within the vcftools module
class vcftools_module_tests (unittest.TestCase):
    # Check vcftools input argument assignment function
    def test_vcf_argument_parser (self):
        # Assign the input argument from the function
        input_arg = vcftools.assign_vcftools_input_arg('example/input/merged_chr1_10000.vcf.gz')
        # Confirm the correct argument has been assigned
        self.assertEqual(input_arg, ['--gzvcf', 'example/input/merged_chr1_10000.vcf.gz'])

    # Check vcftools log output creation
    def test_produce_vcftools_log (self):
        # Create the log output using the function
        vcftools.produce_vcftools_log('Log Test:\n1\n2\n3\n', 'out.logTest')
        # Confirm the correct log output has been created
        self.assertTrue(file_comp('out.logTest.log', 'example/out.logTest.log'))
        # Remove log output file
        self.addCleanup(os.remove, 'out.logTest.log')

    # Check vcftools log for errors
    def test_check_vcftools_for_errors (self):
        # Check the outcome of an output file without errors
        self.assertIsNone(vcftools.check_vcftools_for_errors('Log Test:\n1\n2\n3\nRun Time'))
        # Disable logging module for the following test
        logging.disable(logging.CRITICAL)
        # Check the outcome of an output file with errors
        with self.assertRaises(Exception) as cm:
            vcftools.check_vcftools_for_errors('Log Test:\n1\n2\n3\nError: No Input')
        self.assertEqual(str(cm.exception), 'Error: No Input')

# Run tests for the vcf calc function
class vcf_calc_tests (unittest.TestCase):

    # Check that Fst window function is operating correctly
    def test_Fst_window (self):
        # Run the function with the following arguments
        vcf_calc.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                      '--calc-statistic', 'windowed-weir-fst',
                      '--statistic-window-size', '10000',
                      '--statistic-window-step', '20000',
                      '--model-file', 'example/input/input.model',
                      '--model', '2Pop',
                      '--out-prefix', 'out',
                      '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('out.windowed.weir.fst',
                                  'example/merged_chr1_10000.windowed.weir.fst'))

        # Remove the ouput and log files created by the function
        self.addCleanup(os.remove, 'out.windowed.weir.fst')
        self.addCleanup(os.remove, 'out.windowed.weir.fst.log')

    # Check that the Tajima's D function is operating correctly
    def test_tajimasD (self):
        # Run the function with the following arguments
        vcf_calc.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                      '--calc-statistic', 'TajimaD',
                      '--statistic-window-size', '10000',
                      '--out-prefix', 'out',
                      '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('out.Tajima.D',
                                  'example/merged_chr1_10000.Tajima.D'))

        # Remove the ouput and log files created by the function
        self.addCleanup(os.remove, 'out.Tajima.D')
        self.addCleanup(os.remove, 'out.Tajima.D.log')

    # Check that the window nucleotide diversity function is operating correctly
    def test_window_pi (self):
        # Run the function with the following arguments
        vcf_calc.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                      '--calc-statistic', 'window-pi',
                      '--statistic-window-size', '10000',
                      '--statistic-window-step', '20000',
                      '--out-prefix', 'out',
                      '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('out.windowed.pi',
                                  'example/merged_chr1_10000.windowed.pi'))

        # Remove the ouput and log files created by the function
        self.addCleanup(os.remove, 'out.windowed.pi')
        self.addCleanup(os.remove, 'out.windowed.pi.log')

    # Check that the window nucleotide diversity function is operating correctly
    def test_site_pi (self):
        # Run the function with the following arguments
        vcf_calc.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                      '--calc-statistic', 'site-pi',
                      '--out-prefix', 'out',
                      '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('out.sites.pi',
                                  'example/merged_chr1_10000.sites.pi'))

        # Remove the ouput and log files created by the function
        self.addCleanup(os.remove, 'out.sites.pi')
        self.addCleanup(os.remove, 'out.sites.pi.log')

    # Check that the allele frequency function is operating correctly
    def test_freq (self):
        # Run the function with the following arguments
        vcf_calc.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                      '--calc-statistic', 'freq',
                      '--out-prefix', 'out',
                      '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('out.frq', 'example/merged_chr1_10000.frq'))

        # Remove the ouput and log files created by the function
        self.addCleanup(os.remove, 'out.frq')
        self.addCleanup(os.remove, 'out.frq.log')

    # Check that the heterozygosity function is operating correctly
    def test_het_fit (self):
        # Run the function with the following arguments
        vcf_calc.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                      '--calc-statistic', 'het-fit',
                      '--out-prefix', 'out',
                      '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('out.het', 'example/merged_chr1_10000.het.fit'))

        # Remove the ouput and log files created by the function
        self.addCleanup(os.remove, 'out.het')
        self.addCleanup(os.remove, 'out.het.log')

    # Check that the heterozygosity function is operating correctly
    def test_het_fis (self):
        # Run the function with the following arguments
        vcf_calc.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                      '--model-file', 'example/input/input.model',
                      '--model', '2Pop',
                      '--calc-statistic', 'het-fis',
                      '--out-prefix', 'out',
                      '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('out.het', 'example/merged_chr1_10000.het.fis'))

        # Remove the ouput and log files created by the function
        self.addCleanup(os.remove, 'out.het')
        self.addCleanup(os.remove, 'out.het.log')

# Run tests for the vcf sampler function
class vcf_sampler_tests (unittest.TestCase):

    # Confirm the random sampler is operating correctly
    def test_random_sampler (self):
        # Expected output (based on seed value)
        expected_sample = [967, 713, 222, 321, 898,
                           795, 612, 393, 995, 639,
                           572, 650, 289, 167, 127,
                           809, 831, 820, 210, 725,
                           776, 237, 812, 248, 822]

        # Assign random seed
        np.random.seed(1000)

        # Confirm that the output is what is expected
        self.assertEqual(vcf_sampler.random_vcftools_sampler(range(0, 1000), 25), expected_sample)


    # Confirm the uniform sampler is operating correctly
    def test_uniform_sampler (self):
        # Expected output (based on seed value)
        expected_sample = [118,  70, 163, 161,  73,
                           285, 204, 392, 255, 209,
                           480, 513, 419, 443, 593,
                           731, 620, 627, 605, 723,
                           875, 900, 907, 867, 874]

        # Assign random seed
        np.random.seed(1000)

        # Confirm that the output is what is expected
        self.assertEqual(vcf_sampler.uniform_vcftools_sampler(range(0, 1000), 5, 25), expected_sample)

    # Confirm column assignment is operating correctly
    def test_column_assignment (self):
        # Confirm that the output is what is expected (without errors)
        self.assertEqual(vcf_sampler.assign_position_columns(['CHROM', 'NULL', 'BIN_START', 'BIN_END']), (0, 2, 3))

        # Disable logging module for the following test
        logging.disable(logging.CRITICAL)

        # Confirm that the output is what is expected (with errors)
        with self.assertRaises(ValueError) as cm:
            vcf_sampler.assign_position_columns(['CHROM', 'NULL', 'NULL', 'BIN_END'])
        self.assertEqual(str(cm.exception), 'Cannot find BIN_START column(s) in file specified by --statistic-file.')
        # Confirm that the output is what is expected (with errors)
        with self.assertRaises(ValueError) as cm:
            vcf_sampler.assign_position_columns(['CHROM', 'NULL', 'BIN_START', 'NULL'])
        self.assertEqual(str(cm.exception), 'Cannot find BIN_END column(s) in file specified by --statistic-file.')
        # Confirm that the output is what is expected (with errors)
        with self.assertRaises(ValueError) as cm:
            vcf_sampler.assign_position_columns(['NULL', 'NULL', 'BIN_START', 'BIN_END'])
        self.assertEqual(str(cm.exception), 'Cannot find CHROM column(s) in file specified by --statistic-file.')
        # Confirm that the output is what is expected (with errors)
        with self.assertRaises(ValueError) as cm:
            vcf_sampler.assign_position_columns(['CHROM', 'NULL', 'NULL', 'NULL'])
        self.assertEqual(str(cm.exception), 'Cannot find BIN_START, BIN_END column(s) in file specified by --statistic-file.')
        # Confirm that the output is what is expected (with errors)
        with self.assertRaises(ValueError) as cm:
            vcf_sampler.assign_position_columns(['NULL', 'NULL', 'NULL', 'NULL'])
        self.assertEqual(str(cm.exception), 'Cannot find CHROM, BIN_START, BIN_END column(s) in file specified by --statistic-file.')

    # Check that the entire sampler function (using random sampler) is operating correctly
    def test_sampler (self):
        # Run the function with the following arguments
        vcf_sampler.run(['--vcf', 'example/input/merged_chr1_10000.vcf.gz',
                         '--statistic-file', 'example/merged_chr1_10000.windowed.weir.fst',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--overwrite'])

        # Confirm that the output is what is expected
        self.assertTrue(file_comp('sampled_data.tsv', 'example/sampled_data.tsv'))

        # Remove the ouput files created by the function
        self.addCleanup(os.remove, 'sampled_data.tsv')
        self.addCleanup(shutil.rmtree, 'Sample_Files')

# Run tests for admixr
class admixr_tests (unittest.TestCase):

  # Confirm python version
  def test_version (self):

    print(sys.version)

  # Confirm rpy2 is working
  def test_rpy2 (self):

    # Import rpy2
    from rpy2.robjects import r

    # Obtain the mean from R
    r_mean = r('x <- mean(c(100, 75, 50, 25, 0))')[0]

    # Confirm that the output is what is expected
    self.assertEqual(r_mean, 50.0)



if __name__ == "__main__":
    unittest.main(verbosity = 2)
