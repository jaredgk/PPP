import unittest, filecmp, sys, os
import vcftools_wrapper

def compare_to_expected(test_output, expected_output):
    return filecmp.cmp(test_output, expected_output)

class vcftools_tests (unittest.TestCase):
       
    def test_vcf_argument_parser (self):
        vcf_args = vcftools_wrapper.vcf_argument_parser(['--gzvcf', 'example/locus8.vcf.gz'])
        self.assertEqual(vcf_args.gzvcf, ['--gzvcf', 'example/locus8.vcf.gz'])
    
    def test_produce_vcftools_log (self):
        vcftools_wrapper.produce_vcftools_log('Log Test:\n1\n2\n3\n', 'out', '.logTest')
        self.assertTrue(compare_to_expected('out.logTest.log', 'example/locus8.logTest.log'))
        self.addCleanup(os.remove, 'out.logTest.log')
        
    def test_check_vcftools_for_errors (self):
        self.assertTrue(vcftools_wrapper.check_vcftools_for_errors('Run Time'))
        with self.assertRaises(SystemExit) as cm:
            vcftools_wrapper.check_vcftools_for_errors('SegFault')
        self.assertEqual(cm.exception.code, 'Error with vcftools')

    def test_Fst (self):
        vcftools_wrapper.run(['--gzvcf', 'example/locus8.vcf.gz', '--weir-fst-pop', 'example/Paniscus.txt', '--weir-fst-pop', 'example/Troglodytes.txt', '--out', 'out'])
        self.assertTrue(compare_to_expected('out.weir.fst', 'example/locus8.weir.fst'))
        self.addCleanup(os.remove, 'out.weir.fst')
        self.addCleanup(os.remove, 'out.weir.fst.log')
    
    def test_Fst_window (self):
        vcftools_wrapper.run(['--gzvcf', 'example/locus8.vcf.gz', '--weir-fst-pop', 'example/Paniscus.txt', '--weir-fst-pop', 'example/Troglodytes.txt', '--out', 'out', '--fst-window-size', '--fst-window-step'])
        self.assertTrue(compare_to_expected('out.windowed.weir.fst', 'example/locus8.windowed.weir.fst'))
        self.addCleanup(os.remove, 'out.windowed.weir.fst')
        self.addCleanup(os.remove, 'out.windowed.weir.fst.log')

    def test_tajimasD (self): 
        vcftools_wrapper.run(['--gzvcf', 'example/locus8.vcf.gz', '--out', 'out', '--TajimaD'])
        self.assertTrue(compare_to_expected('out.Tajima.D', 'example/locus8.Tajima.D'))
        self.addCleanup(os.remove, 'out.Tajima.D')
        self.addCleanup(os.remove, 'out.Tajima.D.log')
        
    def test_window_pi (self): 
        vcftools_wrapper.run(['--gzvcf', 'example/locus8.vcf.gz', '--out', 'out', '--window-pi'])
        self.assertTrue(compare_to_expected('out.windowed.pi', 'example/locus8.windowed.pi'))
        self.addCleanup(os.remove, 'out.windowed.pi')
        self.addCleanup(os.remove, 'out.windowed.pi.log')
        
    def test_freq (self):    
        vcftools_wrapper.run(['--gzvcf', 'example/locus8.vcf.gz', '--out', 'out', '--freq'])
        self.assertTrue(compare_to_expected('out.frq', 'example/locus8.frq'))
        self.addCleanup(os.remove, 'out.frq')
        self.addCleanup(os.remove, 'out.frq.log')
    
    def test_het (self):    
        vcftools_wrapper.run(['--gzvcf', 'example/locus8.vcf.gz', '--out', 'out', '--het'])
        self.assertTrue(compare_to_expected('out.het', 'example/locus8.het'))
        self.addCleanup(os.remove, 'out.het')
        self.addCleanup(os.remove, 'out.het.log')

if __name__ == "__main__":
    unittest.main()