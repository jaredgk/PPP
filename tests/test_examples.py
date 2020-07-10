import os
import sys
import unittest
import shutil
import tempfile
import logging
import tracemalloc

from test_functions import fileComp, vcfFileComp

import pgpipe.model_creator as model_creator
import pgpipe.vcf_filter as vcf_filter
import pgpipe.vcf_calc as vcf_calc
import pgpipe.stat_sampler as stat_sampler
import pgpipe.vcf_split as vcf_split

tracemalloc.start()

# Run tests for database.py
class test_examples (unittest.TestCase):

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

    # Test model_creator
    def test_model_creator_example_1 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'model_creator_example_1.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'model_creator.example_1.model')

        # Run the command
        model_creator.run(['--model', '1Pop', '--model-pop', '1Pop', 'Paniscus', '--pop-ind', 'Paniscus', 'Pan_paniscus-9731_LB502', '--out', test_output])

        # Check output
        self.assertTrue(fileComp(test_output, exp_output))

    # Test model_creator
    def test_model_creator_example_2 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'model_creator_example_2.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'model_creator.example_2.model')

        # Assign the input files
        pop_file = os.path.join(self.input_path, '2Pops.txt')
        pan_file = os.path.join(self.input_path, 'Paniscus.txt')
        trog_file = os.path.join(self.input_path, 'Troglodytes.txt')

        # Run the command
        model_creator.run(['--model', '2Pop', '--model-pop-file', '2Pop', pop_file, '--pop-ind-file', 'Paniscus', pan_file, '--pop-ind-file', 'Troglodytes', trog_file, '--out', test_output])

        # Check output
        self.assertTrue(fileComp(test_output, exp_output))

    # Test model_creator
    def test_model_creator_example_3 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'model_creator_example_3.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'model_creator.example_3.model')

        # Assign the input files
        model_file = os.path.join(self.input_path, 'input.model')

        # Run the command
        model_creator.run(['--model-file', model_file, '--update-model', '2Pop', '--model-pop', '2Pop', 'Schweinfurthii', '--model-rm-pop', '2Pop', 'Troglodytes', '--out', test_output])

        # Check output
        self.assertTrue(fileComp(test_output, exp_output))

    # Test model_creator
    def test_model_creator_example_4 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'model_creator_example_4.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'model_creator.example_4.model')

        # Assign the input files
        model_file = os.path.join(self.input_path, 'input.model')

        # Run the command
        model_creator.run(['--model-file', model_file, '--update-pop', 'Paniscus', '--pop-ind', 'Paniscus', 'Pan_paniscus-Unknown', '--pop-rm-ind', 'Paniscus', 'Pan_paniscus-9731_LB502', '--out', test_output])

        # Check output
        self.assertTrue(fileComp(test_output, exp_output))

    # Test vcf_filter
    def test_vcf_filter_example_1 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'vcf_filter_example_1.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.biallelic.vcf.gz')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Run the function with the following arguments
        vcf_filter.run(['--vcf', vcf_file, '--filter-only-biallelic', '--out-format', 'bcf', '--out', test_output])

        # Confirm that the output is what is expected
        self.assertTrue(vcfFileComp(test_output, 'bcf', exp_output, self.test_dir))

    # Test vcf_filter
    def test_vcf_filter_example_2 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'vcf_filter_example_2.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.include_pos.vcf.gz')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.bcf')

        # Run the function with the following arguments
        vcf_filter.run(['--vcf', vcf_file, '--filter-include-pos', 'chr1:1-1509546', '--out-format', 'vcf.gz', '--out', test_output])

        # Confirm that the output is what is expected
        self.assertTrue(vcfFileComp(test_output, 'vcf.gz', exp_output, self.test_dir))

    # Test vcf_filter
    def test_vcf_filter_example_3 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'vcf_filter_example_3.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.no_indels.vcf.gz')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.indels.vcf.gz')

        # Run the function with the following arguments
        vcf_filter.run(['--vcf', vcf_file, '--filter-exclude-indels', '--out-format', 'vcf.gz', '--out', test_output])

        # Confirm that the output is what is expected
        self.assertTrue(vcfFileComp(test_output, 'vcf.gz', exp_output, self.test_dir))

    # Test vcf_calc
    def test_vcf_calc_example_1 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'vcf_calc_example_1.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.Tajima.D')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz') 

        # Run the function with the following arguments
        vcf_calc.run(['--vcf', vcf_file, '--calc-statistic', 'TajimaD', '--statistic-window-size', '10000', '--out', test_output])

        # Confirm that the output is what is expected
        self.assertTrue(fileComp(test_output, exp_output))

    # Test vcf_calc
    def test_vcf_calc_example_2 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'vcf_calc_example_2.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'merged_chr1_10000.windowed.weir.fst')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Assign the input files
        model_file = os.path.join(self.input_path, 'input.model')

        # Run the function with the following arguments
        vcf_calc.run(['--vcf', vcf_file, '--model-file', model_file, '--model', '2Pop', '--calc-statistic', 'windowed-weir-fst', '--statistic-window-size', '10000', '--statistic-window-step', '10000', '--out', test_output])

        #--vcf input.vcf.gz --model-file input.model --model 2Pop --calc-statistic windowed-weir-fst --statistic-window-size 10000 --statistic-window-step 10000

        # Confirm that the output is what is expected
        self.assertTrue(fileComp(test_output, exp_output))

    # Test stat_sampler
    def test_stat_sampler_example_1 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_sampler_example_1.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'sampled.windowed.weir.fst.tsv')

        # Assign the input stat file
        stat_file = os.path.join(self.input_path, 'merged_chr1_10000.windowed.weir.fst')

        # Run the function with the following arguments
        stat_sampler.run(['--statistic-file', stat_file, '--calc-statistic', 'windowed-weir-fst', '--sampling-scheme', 'random', '--sample-size', '20','--random-seed', '1000', '--out', test_output])

        # Confirm that the output is what is expected
        self.assertTrue(fileComp(test_output, exp_output))

    # Test stat_sampler
    def test_stat_sampler_example_2 (self):

        # Assign the test output file
        test_output = os.path.join(self.test_dir, 'test_sampler_example_2.test')

        # Assign the expected output file
        exp_output = os.path.join(self.expected_path, 'sampled.windowed.pi.tsv')

        # Assign the input stat file
        stat_file = os.path.join(self.input_path, 'merged_chr1_10000.windowed.pi')

        # Run the function with the following arguments
        stat_sampler.run(['--statistic-file', stat_file, '--calc-statistic', 'window-pi', '--sampling-scheme', 'uniform', '--uniform-bins', '4', '--sample-size', '20','--random-seed', '1000', '--out', test_output])

        # Confirm that the output is what is expected
        self.assertTrue(fileComp(test_output, exp_output))

    # Test vcf_split
    def test_vcf_split_example_1 (self):

        # Assign the test output file
        test_output_dir = os.path.join(self.test_dir, 'test_vcf_split_example_1')

        # Assign the expected output file
        exp_output_dir = os.path.join(self.expected_path, 'merged_chr1_10000_vcf_split')

        # Assign the input stat file
        stat_file = os.path.join(self.input_path, 'sampled.windowed.weir.fst.tsv')

        # Assign the input vcf
        vcf_file = os.path.join(self.input_path, 'merged_chr1_10000.vcf.gz')

        # Assign the input files
        model_file = os.path.join(self.input_path, 'input.model')

        # Run the function with the following arguments
        vcf_split.run(['--vcf', vcf_file, '--split-file', stat_file, '--split-method', 'statistic-file', '--model-file', model_file, '--model', '2Pop', '--out-dir', test_output_dir])

        # Compare the files within the split dirs
        for file_pos in range(len(os.listdir(test_output_dir))):

            # Assign the file paths
            test_output = os.path.join(test_output_dir, 'out_%s.vcf.gz' % file_pos)
            exp_output = os.path.join(exp_output_dir, 'out_%s.vcf.gz' % file_pos)

            # Confirm that the output is what is expected
            self.assertTrue(vcfFileComp(test_output, 'vcf.gz', exp_output, self.test_dir))





        

if __name__ == "__main__":
    unittest.main(verbosity = 2)