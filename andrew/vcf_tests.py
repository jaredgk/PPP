import unittest, filecmp, sys, os
import vcf_scripts

def compare_to_expected(test_output, expected_output):
    return filecmp.cmp(test_output, expected_output)

class vcftools_tests (unittest.TestCase):
       
    def test_vcf_argument_parser (self):
        vcf_args = vcf_scripts.vcf_argument_parser(['--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz'])
        self.assertEqual(vcf_args.input, ['--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz'])
    
    def test_produce_vcftools_log (self):
        vcf_scripts.produce_vcftools_log('Log Test:\n1\n2\n3\n', 'out', '.logTest')
        self.assertTrue(compare_to_expected('andrew' + '/' + 'out.logTest.log', 'andrew' + '/' + 'example/locus8.logTest.log'))
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.logTest.log')
    
    def test_return_basic_args (self):
        vcf_args = vcf_scripts.vcf_argument_parser(['--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz'])
        self.assertEqual(vcf_scripts.return_basic_args(vcf_args), ['vcftools', '--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz', '--out', 'out'])
    
    def test_check_vcftools_for_errors (self):
        self.assertTrue(vcf_scripts.check_vcftools_for_errors('Run Time'))
        with self.assertRaises(SystemExit) as cm:
            vcf_scripts.check_vcftools_for_errors('SegFault')
        self.assertEqual(cm.exception.code, 'Error with vcftools')

    def test_calculate_Fst (self):
        vcf_scripts.calculate_Fst(['--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz', '--weir-fst-pops', 'andrew' + '/' + 'example/Paniscus.txt', 'andrew' + '/' + 'example/Troglodytes.txt', '--out', 'out'])
        self.assertTrue(compare_to_expected('andrew' + '/' + 'out.windowed.weir.fst', 'andrew' + '/' + 'example/locus8.windowed.weir.fst'))
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.windowed.weir.fst')
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.windowed.weir.fst.log')
    
    def test_calculate_tajimasD (self): 
        vcf_scripts.calculate_tajimasD(['--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz', '--out', 'out'])
        self.assertTrue(compare_to_expected('out.Tajima.D', 'andrew' + '/' + 'example/locus8.Tajima.D'))
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.Tajima.D')
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.Tajima.D.log')
    
    def test_calculate_pi (self): 
        vcf_scripts.calculate_pi(['--gzvcf', 'example/locus8.vcf.gz', '--out', 'out'])
        self.assertTrue(compare_to_expected('andrew' + '/' + 'out.windowed.pi', 'andrew' + '/' + 'example/locus8.windowed.pi'))
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.windowed.pi')
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.windowed.pi.log')
    
    def test_calculate_af(self):    
        vcf_scripts.calculate_af(['--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz', '--out', 'out'])
        self.assertTrue(compare_to_expected('andrew' + '/' + 'out.frq', 'andrew' + '/' + 'example/locus8.frq'))
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.frq')
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.frq.log')
    
    def test_calculate_het(self):    
        vcf_scripts.calculate_heterozygosity(['--gzvcf', 'andrew' + '/' + 'example/locus8.vcf.gz', '--out', 'out'])
        self.assertTrue(compare_to_expected('andrew' + '/' + 'out.het', 'andrew' + '/' + 'example/locus8.het'))
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.het')
        self.addCleanup(os.remove, 'andrew' + '/' + 'out.het.log')

if __name__ == "__main__":
    unittest.main()