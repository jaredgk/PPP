import unittest
import pysam
import filecmp
import os
from vcf_ref_to_seq import getMaxAlleleLength, getFastaFilename, \
    vcf_to_seq, createParser, validateFiles
from gene_region import RegionList, Region


class funcTest(unittest.TestCase):

    def test_getMaxAlleleLength(self):
        l = ('AAA', 'A', 'AAAAAAA')
        self.assertEqual(getMaxAlleleLength(l), 7)

    def test_getFastaFilename(self):
        fn = 'test.vcf.gz'
        self.assertEqual(getFastaFilename(fn), ('test.fasta', '.vcf.gz'))

    def testPathValidator(self):
        parser = createParser()
        args = ['--vcf', 'example/doesntexist.vcf.gz',
                '--ref', 'example/human_g1k_chr11.fasta',
                '--gr', 'example/snp_region.txt']
        pa = parser.parse_args(args)
        self.assertRaises(ValueError, validateFiles, pa)


class geneRegionTest(unittest.TestCase):

    def test_RL_collist(self):
        collist = [2, 4, 6]
        rl = RegionList('example/gr_ex_multicolumn.txt', collist=collist)
        r = rl.regions[0]
        set_r = Region(99999, 100099, '11')
        self.assertEqual(r, set_r)


class snpTest(unittest.TestCase):

    def test_generateSequence_snp(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--gr', 'example/snp_region.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def test_generateSequence_insert(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--gr', 'example/insert_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.insert.example.fasta'), True)

    def test_generateSequence_del(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--gr', 'example/del_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.del.example.fasta'), True)

    def test_generateSequence_multi(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--gr', 'example/multi_inregion.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.multi.example.fasta'), True)

    def test_subsample(self):
        vcf_to_seq(['--vcf', 'example/chr11.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--gr', 'example/snp_region.txt',
             '--subsamp_list', 'example/subsample_list.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def tearDown(self):
        try:
            os.remove('example/chr11.subsamples.fasta')
            #figure out how to do this if test fails
        except:
            pass
        try:
            os.remove('example/chr11.fasta')
        except:
            pass


if __name__ == "__main__":
    unittest.main()
