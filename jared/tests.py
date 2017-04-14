import unittest
import pysam
import parser
import filecmp
import os
from vcf_ref_to_seq import getMaxAlleleLength, getFastaFilename, \
    vcf_to_seq, createParser, validateFiles
from vcf_from_regions import vcf_region_write
from gene_region import RegionList, Region

def tryRemove(filename):
    try:
        os.remove(filename)
    except:
        pass

class funcTest(unittest.TestCase):

    def test_getMaxAlleleLength(self):
        l = ('AAA', 'A', 'AAAAAAA')
        self.assertEqual(getMaxAlleleLength(l), 7)

    def test_getFastaFilename(self):
        fn = 'test.vcf.gz'
        p = createParser()
        args = p.parse_args(['--vcf', 'test.vcf.gz'])
        self.assertEqual(getFastaFilename(args), ('test.fasta', 'vcf.gz'))

    def testPathValidator(self):
        parser = createParser()
        args = ['--vcf', 'example/doesntexist.vcf.gz',
                '--ref', 'example/human_g1k_chr11.fasta',
                '--rl', 'example/snp_region.txt']
        pa = parser.parse_args(args)
        self.assertRaises(ValueError, validateFiles, pa)


class geneRegionTest(unittest.TestCase):

    def test_RL_collist(self):
        collist = '2,4,6'
        rl = RegionList('example/gr_ex_multicolumn.txt', colstr=collist)
        r = rl.regions[0]
        set_r = Region(100000, 100100, '11')
        self.assertEqual(r, set_r)

    def test_RL_sort(self):
        rl = RegionList('example/unsorted_regions.txt')
        sl = RegionList('example/sorted_regions.txt')
        rl.regions.sort()
        self.assertTrue(rl.regions == sl.regions)


class snpTest(unittest.TestCase):

    def test_generateSequence_snp(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--rl', 'example/snp_region.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def test_generateSequence_insert(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--rl', 'example/insert_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.insert.example.fasta'), True)

    def test_generateSequence_del(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--rl', 'example/del_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.del.example.fasta'), True)

    def test_generateSequence_multi(self):
        vcf_to_seq(['--vcf', 'example/chr11.subsamples.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--rl', 'example/multi_inregion.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.multi.example.fasta'), True)

    def test_subsample(self):
        vcf_to_seq(['--vcf', 'example/chr11.vcf.gz',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--rl', 'example/snp_region.txt',
             '--subsamp_list', 'example/subsample_list.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def test_generateSequence_unzipped(self):
        vcf_to_seq(['--vcf', 'example/chr11.unzipped.vcf',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--rl', 'example/snp_region.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.unzipped.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def test_generateSequence_uz_adj(self):
        vcf_to_seq(['--vcf', 'example/chr11.unzipped.vcf',
             '--ref', 'example/human_g1k_chr11.fasta',
             '--rl', 'example/adj_regions.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.unzipped.fasta',
                         'example/chr11.adjacent.fasta'), True)

    def test_generateSequence_selfcompress(self):
        tryRemove('example/chr11.unzipped.vcf.gz')
        vcf_to_seq(['--vcf', 'example/chr11.unzipped.vcf',
                '--ref', 'example/human_g1k_chr11.fasta',
                '--rl', 'example/snp_region.txt', '--compress-vcf'])
        self.assertTrue(filecmp.cmp('example/chr11.unzipped.fasta',
                        'example/chr11.snp.example.fasta'))

    def tearDown(self):
        tryRemove('example/chr11.subsamples.fasta')
        tryRemove('example/chr11.fasta')
        tryRemove('example/chr11.unzipped.fasta')

class reduceTest(unittest.TestCase):

    def test_vcf_region_write_simple(self):
        vcf_region_write(['example/chr11.subsamples.vcf.gz',
                         '--rl','example/snp_region.txt',
                         '--output','example/chr11.temp.vcf.gz'])
        self.assertTrue(filecmp.cmp('example/chr11.snp.sr.vcf.gz',
                        'example/chr11.temp.vcf.gz'))

    def tearDown(self):
        tryRemove('example/chr11.temp.vcf.gz')



if __name__ == "__main__":
    unittest.main()
