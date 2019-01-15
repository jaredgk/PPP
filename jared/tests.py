import unittest
import pysam
import parser
import filecmp
import os
import sys
from vcf_ref_to_seq import getMaxAlleleLength, getFastaFilename, \
    vcf_to_seq, createParser, validateFiles
from vcf_from_regions import vcf_region_write
from vcf_reader_func import checkFormat, VcfReader
from gene_region import RegionList, Region
from logging_module import initLogger
from four_gamete_pysam import sample_fourgametetest_intervals
from vcf_to_ima import vcf_to_ima
from find_intergenic_bed import get_intergenic
from get_nonmissing_chunks import regionsWithData

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
        args = p.parse_args(['test.vcf.gz','test.fasta','test.reg'])
        self.assertEqual(getFastaFilename(args), 'test.fasta')

    def testPathValidator(self):
        parser = createParser()
        args = ['example/doesntexist.vcf.gz',
                'example/human_g1k_chr11.fasta',
                'example/snp_region.txt']
        pa = parser.parse_args(args)
        self.assertRaises(ValueError, validateFiles, pa)

class intergenicTest(unittest.TestCase):
    def test_oneidx(self):
        args = ['--bed','example/fib_in.txt',
                '--out','example/fib_t.txt']
        get_intergenic(args)
        self.assertTrue(filecmp.cmp('example/fib_t.txt','example/fib_oneidx.txt'))
    def test_zeroho(self):
        args = ['--bed','example/fib_in.txt','--out','example/fib_t.txt','--zero-ho']
        get_intergenic(args)
        self.assertTrue(filecmp.cmp('example/fib_t.txt','example/fib_zeroho.txt'))

class missingTest(unittest.TestCase):
    def test_default(self):
        args = ['--vcf','example/macaca_missingdata.vcf.gz',
                '--out','example/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('example/gmd_test.txt','example/gmd_test_default.txt'))

    def test_extend(self):
        args = ['--vcf','example/macaca_missingdata.vcf.gz',
                '--extend-regions',
                '--out','example/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('example/gmd_test.txt','example/gmd_test_extend.txt'))

    def test_allowOneMiss(self):
        args = ['--vcf','example/macaca_missingdata.vcf.gz',
                '--missing-count','1',
                '--out','example/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('example/gmd_test.txt','example/gmd_test_miss1.txt'))

    def test_allsize(self):
        args = ['--vcf','example/macaca_missingdata.vcf.gz',
                '--size','0',
                '--out','example/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('example/gmd_test.txt','example/gmd_test_size0.txt'))


class geneRegionTest(unittest.TestCase):

    def getARecord(self):
        vcfr = VcfReader('example/chr11.subsamples.vcf.gz')
        return vcfr.prev_last_rec

    def test_RL_collist(self):
        collist = '2,4,6'
        rl = RegionList('example/gr_ex_multicolumn.txt', colstr=collist, zeroho=True)
        r = rl.regions[0]
        set_r = Region(100000, 100100, '11')
        self.assertEqual(r, set_r)

    def test_RL_sort(self):
        rl = RegionList('example/unsorted_regions.txt')
        sl = RegionList('example/sorted_regions.txt')
        rl.regions.sort()
        self.assertTrue(rl.regions == sl.regions)

    def test_RL_charsort(self):
        rl = RegionList('example/unsorted_regions.txt',sortmethod="string")
        sl = RegionList('example/charsorted_regions.txt')
        rl.regions.sort()
        self.assertTrue(rl.regions == sl.regions)


    def test_RL_overlap(self):
        rl = RegionList('example/overlap_regions.txt')
        self.assertTrue(rl.hasOverlap())

    def test_CR_in(self):
        rec = self.getARecord()
        rl = RegionList(genestr="11:142500:142600")
        self.assertTrue(rl.regions[0].containsRecord(rec) == 'in')

    def test_CR_before(self):
        rec = self.getARecord()
        rl = RegionList(genestr="11:142531:142600")
        val = rl.regions[0].containsRecord(rec)
        sys.stderr.write(str(val)+'\n')
        self.assertTrue(rl.regions[0].containsRecord(rec) == 'before')

    def test_CR_after(self):
        rec = self.getARecord()
        rl = RegionList(genestr="11:142500:142529")
        val = rl.regions[0].containsRecord(rec)
        self.assertTrue(val == 'after')


class snpTest(unittest.TestCase):

    def test_generateSequence_snp(self):
        vcf_to_seq(['example/chr11.subsamples.vcf.gz',
             'example/human_g1k_chr11.fasta',
             'example/snp_region.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def test_generateSequence_insert(self):
        vcf_to_seq(['example/chr11.subsamples.vcf.gz',
             'example/human_g1k_chr11.fasta',
             'example/insert_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.insert.example.fasta'), True)

    def test_generateSequence_del(self):
        vcf_to_seq(['example/chr11.subsamples.vcf.gz',
             'example/human_g1k_chr11.fasta',
             'example/del_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.del.example.fasta'), True)

    def test_generateSequence_multi(self):
        vcf_to_seq(['example/chr11.subsamples.vcf.gz',
             'example/human_g1k_chr11.fasta',
             'example/multi_inregion.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.subsamples.fasta',
                         'example/chr11.multi.example.fasta'), True)

    def test_subsample(self):
        vcf_to_seq(['example/chr11.vcf.gz',
             'example/human_g1k_chr11.fasta',
             'example/snp_region.txt',
             '--subsamp-list', 'example/subsample_list.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def test_generateSequence_unzipped(self):
        vcf_to_seq(['example/chr11.unzipped.vcf',
             'example/human_g1k_chr11.fasta',
             'example/snp_region.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.unzipped.fasta',
                         'example/chr11.snp.example.fasta'), True)

    def test_generateSequence_uz_adj(self):
        vcf_to_seq(['example/chr11.unzipped.vcf',
             'example/human_g1k_chr11.fasta',
            'example/adj_regions.txt'])
        self.assertEqual(filecmp.cmp('example/chr11.unzipped.fasta',
                         'example/chr11.adjacent.fasta'), True)

    def test_generateSequence_selfcompress(self):
        tryRemove('example/chr11.unzipped.vcf.gz')
        vcf_to_seq(['example/chr11.unzipped.vcf',
                'example/human_g1k_chr11.fasta',
                'example/snp_region.txt', '--compress-vcf'])
        self.assertTrue(filecmp.cmp('example/chr11.unzipped.fasta',
                        'example/chr11.snp.example.fasta'))

    def test_generateSequence_mismatchref(self):
        args = ['example/chr11.subsamples.vcf.gz',
             'example/human_g1k_messy.fasta',
             'example/snp_region.txt']

        self.assertRaises(Exception, vcf_to_seq, args)

    def tearDown(self):
        tryRemove('example/chr11.subsamples.fasta')
        tryRemove('example/chr11.fasta')
        tryRemove('example/chr11.unzipped.fasta')

class checkCompressionTest(unittest.TestCase):

    def test_checkFormat_nozip(self):
        comp = checkFormat('example/chr11.unzipped.vcf')
        self.assertTrue(comp == 'vcf')

    def test_checkFormat_bgzip(self):
        comp = checkFormat('example/chr11.subsamples.vcf.gz')
        self.assertTrue(comp == 'bgzip')

    def test_checkFormat_gzip(self):
        comp = checkFormat('example/chr11.gzipped.vcf.gz')
        self.assertTrue(comp=='gzip')

    def test_checkFormat_bcf(self):
        comp = checkFormat('example/chr11.test.bcf')
        self.assertTrue(comp=='bcf')

    def test_checkFormat_other(self):
        comp = checkFormat('example/plain.txt')
        self.assertTrue(comp=="other")

    def test_checkFormat_otherzip(self):
        comp = checkFormat('example/plain.txt.gz')
        self.assertTrue(comp=="other")

class reduceTest(unittest.TestCase):

    def test_vcf_region_write_simple(self):
        vcf_region_write(['example/chr11.subsamples.vcf.gz',
                         '--rl','example/snp_region.txt',
                         '--output','example/chr11.test.vcf.gz'])
        self.assertTrue(filecmp.cmp('example/chr11.snp.sr.vcf.gz',
                        'example/chr11.test.vcf.gz'))

    def test_vcf_region_write_unzipped(self):
        vcf_region_write(['example/chr11.unzipped.vcf',
                          '--rl','example/snp_region.txt',
                          '--output','example/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('example/chr11.snp.sr.vcf',
                        'example/chr11.test.vcf'))

    def test_vcf_region_write_first(self):
        vcf_region_write(['example/chr11.unzipped.vcf',
                          '--rl','example/first_region.txt',
                          '--output','example/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('example/chr11.first.vcf',
                        'example/chr11.test.vcf'))

    def test_vcf_region_write_last(self):
        vcf_region_write(['example/chr11.unzipped.vcf',
                          '--rl','example/last_region.txt',
                          '--output','example/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('example/chr11.last.vcf',
                        'example/chr11.test.vcf'))

    def test_vcf_region_write_multi(self):
        vcf_region_write(['example/chr11.unzipped.vcf',
                          '--rl','example/multi_big_region.txt',
                          '--output','example/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('example/chr11.multi.big.vcf',
                        'example/chr11.test.vcf'))

    def test_vcf_region_write_adj(self):
        vcf_region_write(['example/chr11.unzipped.vcf',
                          '--rl','example/adj_big_region.txt',
                          '--output','example/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('example/chr11.adj.big.vcf',
                        'example/chr11.test.vcf'))

    def tearDown(self):
        tryRemove('example/chr11.test.vcf.gz')
        tryRemove('example/chr11.test.vcf')



class configTest(unittest.TestCase):

    def test_vcf_to_seq_config(self):
        vcf_to_seq(['example/chr11.subsamples.vcf.gz',
                    'example/human_g1k_chr11.fasta',
                    'example/insert_region.txt',
                    '--conf','example/config.config'])
        self.assertTrue(filecmp.cmp('example/chr11.insert.example.fasta',
                        'example/chr11.subsamples.fasta'))

    def test_vcf_to_seq_config_required(self):
        vcf_to_seq(['--conf', 'example/withreq.conf'])
        self.assertTrue(filecmp.cmp('example/chr11.snp.example.fasta',
                        'example/chr11.subsamples.fasta'))

class fourgameteTest(unittest.TestCase):

    def test_hk(self):
        region_list = [[[192385, 196945], [198986, 202254], [202547, 202813],
         [206421, 207401], [207400, 215365], [217814, 218641]]]
        test_list = sample_fourgametetest_intervals(["--vcfname","example/chr11subsamples4gtest.vcf.gz","--hk","--retl"])
        self.assertEqual(region_list, test_list)

    def test_comp(self):
        region_list = [[[192385, 196943], [192386, 202252], [198987, 202546],
         [201585, 202811], [202548, 207399], [206422, 215363], [207401, 218140],
          [207463, 218639], [217815, 221584]]]
        test_list = sample_fourgametetest_intervals(["--vcfname","example/chr11subsamples4gtest.vcf.gz","--4gcompat","--retl",'--ovlpi','--ovlps'])
        self.assertEqual(region_list,test_list)

    def test_vcf(self):
        sample_fourgametetest_intervals(["--vcfname","example/chr11.subsamples.vcf.gz",'--4gcompat','--reti','--out','example/chr11.4g.vcf','--numinf','2'])
        self.assertTrue(filecmp.cmp('example/chr11.4gtest1.vcf','example/chr11.4g.vcf'))
        tryRemove('example/chr11.4g.vcf')

    def test_consecutive(self):
        sample_fourgametetest_intervals(["--vcfname","example/chr11.4gsmall.vcf",'--4gcompat','--reti','--out','example/chr11.4g.vcf','--numinf','2','--ovlpi','--ovlps'])
        self.assertTrue(filecmp.cmp('example/chr11.4gtest2.vcf','example/chr11.4g.vcf'))
        tryRemove('example/chr11.4g.vcf')


class imaTest(unittest.TestCase):
    def test_ima(self):
        vcf_to_ima(['--vcf','example/chr11.subsamples.vcf.gz',
                    '--ref','example/human_g1k_chr11.fasta',
                    '--bed','example/snp_region.txt',
                    '--pop','example/testmodel.model',
                    '--zero-ho',
                    '--output','example/chr11.subsamples.ima.u'])
        self.assertTrue(filecmp.cmp('example/chr11.ima.u',
                        'example/chr11.subsamples.ima.u'))



if __name__ == "__main__":
    initLogger(filelevel="ERROR",streamlevel="ERROR")
    unittest.main(verbosity=2)
