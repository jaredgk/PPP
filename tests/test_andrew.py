import unittest
import filecmp
import sys
import os
import logging
import shutil
import tempfile
import numpy as np

# Import scripts to test
import pgpipe.vcftools as vcftools
import pgpipe.vcf_calc as vcf_calc
import pgpipe.stat_sampler as stat_sampler

# Used to compare two files, returns bool
def file_comp(test_output, expected_output):
    return filecmp.cmp(test_output, expected_output)

# Run tests for the functions within the vcftools module
class vcftools_module_tests (unittest.TestCase):

    @classmethod
    def setUpClass (cls):

        # Supress warnings
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)

        # Create a temporary directory
        cls.test_dir = tempfile.mkdtemp()

        # Assign the script directory
        cls.script_dir = os.path.dirname(os.path.realpath(__file__))

        # Assign the expected output directory
        cls.expected_dir = 'exp_output'

        # Assign the expected path
        cls.expected_path = os.path.join(cls.script_dir, cls.expected_dir)

        # Assign the expected path
        cls.input_path = os.path.join(cls.script_dir, 'input')

    @classmethod
    def tearDownClass (cls):

        # Remove the test directory after the tests
        shutil.rmtree(cls.test_dir)

    # Check vcftools input argument assignment function
    def test_vcf_argument_parser (self):

        # Assign the input argument from the function
        input_arg = vcftools.assign_vcftools_input_arg('example/input/merged_chr1_10000.vcf.gz')

        # Confirm the correct argument has been assigned
        self.assertEqual(input_arg, ['--gzvcf', 'example/input/merged_chr1_10000.vcf.gz'])

    # Check vcftools log output creation
    def test_produce_vcftools_log (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_produce_vcftools_log')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'out.logTest.log')

        # Create the log output using the function
        vcftools.produce_vcftools_log('Log Test:\n1\n2\n3\n', test_output)

        # Confirm the correct log output has been created
        self.assertTrue(file_comp(test_output + '.log', exp_output))

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

    @classmethod
    def setUpClass (cls):

        # Supress warnings
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)

        # Create a temporary directory
        cls.test_dir = tempfile.mkdtemp()

        # Assign the script directory
        cls.script_dir = os.path.dirname(os.path.realpath(__file__))

        # Assign the expected output directory
        cls.expected_dir = 'exp_output'

        # Assign the expected path
        cls.expected_path = os.path.join(cls.script_dir, cls.expected_dir)

        # Assign the expected path
        cls.input_path = os.path.join(cls.script_dir, 'input')

    @classmethod
    def tearDownClass (cls):

        # Remove the test directory after the tests
        shutil.rmtree(cls.test_dir)

    # Check that Fst window function is operating correctly
    def test_Fst_window (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_Fst_window.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.windowed.weir.fst')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Assign the input files
        model_file = os.path.join(self.input_path, 'input.model')

        # Run the function with the following arguments
        vcf_calc.run(vcf = vcf_file,
                     calc_statistic = 'windowed-weir-fst',
                     statistic_window_size = 10000,
                     statistic_window_step = 10000,
                     model_file = model_file,
                     model = '2Pop',
                     out = test_output)

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(test_output, exp_output))

    # Check that the Tajima's D function is operating correctly
    def test_tajimasD (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_tajimasD.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.Tajima.D')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Run the function with the following arguments
        vcf_calc.run(vcf = vcf_file,
                     calc_statistic = 'TajimaD',
                     statistic_window_size = 10000,
                     out = test_output)

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(test_output, exp_output))

    # Check that the window nucleotide diversity function is operating correctly
    def test_window_pi (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_window_pi.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.windowed.pi')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Run the function with the following arguments
        vcf_calc.run(vcf = vcf_file,
                     calc_statistic = 'window-pi',
                     statistic_window_size = 10000,
                     statistic_window_step = 20000,
                     out = test_output)

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(test_output, exp_output))

    # Check that the window nucleotide diversity function is operating correctly
    def test_site_pi (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_site_pi.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.sites.pi')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Run the function with the following arguments
        vcf_calc.run(vcf = vcf_file,
                     calc_statistic =  'site-pi',
                     out = test_output)

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(test_output, exp_output))

    # Check that the allele frequency function is operating correctly
    def test_freq (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_freq.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.frq')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Run the function with the following arguments
        vcf_calc.run(vcf = vcf_file,
                     calc_statistic = 'freq',
                     out = test_output)

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(test_output, exp_output))

    # Check that the heterozygosity function is operating correctly
    def test_het_fit (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_het_fit.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.het.fit')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Run the function with the following arguments
        vcf_calc.run(vcf = vcf_file,
                     calc_statistic = 'het-fit',
                     out = test_output)

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(test_output, exp_output))

    # Check that the heterozygosity function is operating correctly
    def test_het_fis (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_het_fis.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.het.fis')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Assign the input files
        model_file = os.path.join(self.input_path, 'input.model')

        # Run the function with the following arguments
        vcf_calc.run(vcf = vcf_file,
                     model_file = model_file,
                     model = '2Pop',
                     calc_statistic = 'het-fis',
                     out = test_output)

        # Confirm that the output is what is expected
        self.assertTrue(file_comp(test_output, exp_output))

# Run tests for the vcf sampler function
class vcf_sampler_tests (unittest.TestCase):

    @classmethod
    def setUpClass (cls):

        # Supress warnings
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)

        # Create a temporary directory
        cls.test_dir = tempfile.mkdtemp()

        # Assign the script directory
        cls.script_dir = os.path.dirname(os.path.realpath(__file__))

        # Assign the expected output directory
        cls.expected_dir = 'exp_output'

        # Assign the expected path
        cls.expected_path = os.path.join(cls.script_dir, cls.expected_dir)

        # Assign the expected path
        cls.input_path = os.path.join(cls.script_dir, 'input')

    @classmethod
    def tearDownClass (cls):

        # Remove the test directory after the tests
        shutil.rmtree(cls.test_dir)

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
        self.assertEqual(stat_sampler.random_vcftools_sampler(range(0, 1000), 25), expected_sample)

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
        self.assertEqual(stat_sampler.uniform_vcftools_sampler(range(0, 1000), 5, 25), expected_sample)

    # Confirm column assignment is operating correctly
    def test_column_assignment (self):

        # Confirm that the output is what is expected (without errors)
        self.assertEqual(stat_sampler.assign_statistic_column(['CHROM', 'NULL', 'BIN_START', 'BIN_END', 'MEAN_FST'], 'windowed-weir-fst'), 'MEAN_FST')
        self.assertEqual(stat_sampler.assign_statistic_column(['CHROM', 'NULL', 'BIN_START', 'BIN_END', 'TajimaD'], 'TajimaD'), 'TajimaD')
        self.assertEqual(stat_sampler.assign_statistic_column(['CHROM', 'NULL', 'BIN_START', 'BIN_END', 'PI'], 'window-pi'), 'PI')

# Run tests for admixr
class admixr_tests (unittest.TestCase):

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
