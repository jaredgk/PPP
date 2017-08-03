import unittest, filecmp, sys, os
import vcftools
import vcf_calc

def compare_to_expected(test_output, expected_output):
    return filecmp.cmp(test_output, expected_output)

class vcftools_tests (unittest.TestCase):

    def test_vcf_argument_parser (self):
        input_arg = vcftools.assign_vcftools_input_arg('example/locus8.vcf.gz')
        self.assertEqual(input_arg, ['--gzvcf', 'example/locus8.vcf.gz'])

    def test_produce_vcftools_log (self):
        vcftools.produce_vcftools_log('Log Test:\n1\n2\n3\n', 'out', 'logTest')
        self.assertTrue(compare_to_expected('out.logTest.log', 'example/locus8.logTest.log'))
        self.addCleanup(os.remove, 'out.logTest.log')

    def test_check_vcftools_for_errors (self):
        self.assertTrue(vcftools.check_vcftools_for_errors('Log Test:\n1\n2\n3\nRun Time'))
        with self.assertRaises(Exception) as cm:
            vcftools.check_vcftools_for_errors('Log Test:\n1\n2\n3\nError: No Input')
        self.assertEqual(cm.exception.code, 'Error: No Input')

    def test_Fst_window (self):
        vcf_calc.run(['example/locus8.vcf.gz', '--calc-statistic', 'windowed-weir-fst', '--pop-file', 'example/Paniscus.txt', '--pop-file', 'example/Troglodytes.txt', '--out', 'out'])
        self.assertTrue(compare_to_expected('out.windowed.weir.fst', 'example/locus8.windowed.weir.fst'))
        self.addCleanup(os.remove, 'out.windowed.weir.fst')
        self.addCleanup(os.remove, 'out.windowed.weir.fst.log')

    def test_tajimasD (self):
        vcf_calc.run(['example/locus8.vcf.gz', '--calc-statistic', 'TajimaD', '--out', 'out'])
        self.assertTrue(compare_to_expected('out.Tajima.D', 'example/locus8.Tajima.D'))
        self.addCleanup(os.remove, 'out.Tajima.D')
        self.addCleanup(os.remove, 'out.Tajima.D.log')

    def test_window_pi (self):
        vcf_calc.run(['example/locus8.vcf.gz', '--calc-statistic', 'pi', '--out', 'out'])
        self.assertTrue(compare_to_expected('out.windowed.pi', 'example/locus8.windowed.pi'))
        self.addCleanup(os.remove, 'out.windowed.pi')
        self.addCleanup(os.remove, 'out.windowed.pi.log')

    def test_freq (self):
        vcf_calc.run(['example/locus8.vcf.gz', '--calc-statistic', 'freq', '--out', 'out'])
        self.assertTrue(compare_to_expected('out.frq', 'example/locus8.frq'))
        self.addCleanup(os.remove, 'out.frq')
        self.addCleanup(os.remove, 'out.frq.log')

    def test_het (self):
        vcf_calc.run(['example/locus8.vcf.gz', '--calc-statistic', 'het', '--out', 'out'])
        self.assertTrue(compare_to_expected('out.het', 'example/locus8.het'))
        self.addCleanup(os.remove, 'out.het')
        self.addCleanup(os.remove, 'out.het.log')

if __name__ == "__main__":
    unittest.main()
