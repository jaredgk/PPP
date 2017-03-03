import unittest
import pysam
import filecmp
import os
from vcf_ref_to_seq import getMaxAlleleLength,getFastaFilename,main



class FirstTestOfGmal(unittest.TestCase):

    def test_getMaxAlleleLength(self):
        l = ('AAA','A','AAAAAAA')
        self.assertEqual(getMaxAlleleLength(l),7)

    def test_getFastaFilename(self):
        fn = 'test.vcf.gz'
        self.assertEqual(getFastaFilename(fn),'test.fasta')


class snpTest(unittest.TestCase):

    def test_generateSequence_snp(self):
        main(['--vcf','example/chr11.subsamples.vcf.gz','--ref',
            'example/human_g1k_chr11.fasta',
            '--gr','example/snp_region.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.snp.example.fasta'),True)

    def test_generateSequence_insert(self):
        main(['--vcf','example/chr11.subsamples.vcf.gz',
            '--ref','example/human_g1k_chr11.fasta',
            '--gr','example/insert_region.txt','--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.insert.example.fasta'),True)

    def test_generateSequence_del(self):
        main(['--vcf','example/chr11.subsamples.vcf.gz',
            '--ref','example/human_g1k_chr11.fasta',
            '--gr','example/del_region.txt','--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                        'example/chr11.del.example.fasta'),True)

    def test_generateSequence_multi(self):
        main(['--vcf','example/chr11.subsamples.vcf.gz',
            '--ref','example/human_g1k_chr11.fasta',
            '--gr','example/multi_inregion.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                        'example/chr11.multi.example.fasta'),True)

    def tearDown(self):
        try:
            os.remove('example/chr11.subsamples.fasta')
            #figure out how to do this if test fails
        except:
            pass




if __name__ == "__main__":
    unittest.main()
