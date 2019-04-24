import unittest
import pysam
import parser
import filecmp
import os
import sys
from pgpipe.vcf_ref_to_seq import getMaxAlleleLength, getFastaFilename, \
    vcf_to_seq, createParser, validateFiles
from pgpipe.vcf_from_regions import vcf_region_write
from pgpipe.vcf_reader_func import checkFormat, VcfReader
from pgpipe.gene_region import RegionList, Region
from pgpipe.logging_module import initLogger
from pgpipe.four_gamete_pysam import sample_fourgametetest_intervals
from pgpipe.vcf_to_ima import vcf_to_ima
from pgpipe.find_intergenic_bed import get_intergenic
from pgpipe.get_nonmissing_chunks import regionsWithData
from pgpipe import vcf_phase

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
        args = ['input/doesntexist.vcf.gz',
                'input/human_g1k_chr11.fasta',
                'input/snp_region.txt']
        pa = parser.parse_args(args)
        self.assertRaises(ValueError, validateFiles, pa)

class intergenicTest(unittest.TestCase):
    def test_oneidx(self):
        args = ['--bed','input/fib_in.txt',
                '--out','input/fib_t.txt']
        get_intergenic(args)
        self.assertTrue(filecmp.cmp('input/fib_t.txt','input/fib_oneidx.txt'))
    def test_zeroho(self):
        args = ['--bed','input/fib_in.txt','--out','input/fib_t.txt','--zero-ho']
        get_intergenic(args)
        self.assertTrue(filecmp.cmp('input/fib_t.txt','input/fib_zeroho.txt'))

class missingTest(unittest.TestCase):
    def test_default(self):
        args = ['--vcf','input/macaca_missingdata.vcf.gz',
                '--out','input/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('input/gmd_test.txt','input/gmd_test_default.txt'))

    def test_extend(self):
        args = ['--vcf','input/macaca_missingdata.vcf.gz',
                '--extend-regions',
                '--out','input/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('input/gmd_test.txt','input/gmd_test_extend.txt'))

    def test_allowOneMiss(self):
        args = ['--vcf','input/macaca_missingdata.vcf.gz',
                '--missing-count','1',
                '--out','input/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('input/gmd_test.txt','input/gmd_test_miss1.txt'))

    def test_allsize(self):
        args = ['--vcf','input/macaca_missingdata.vcf.gz',
                '--size','0',
                '--out','input/gmd_test.txt']
        regionsWithData(args)
        self.assertTrue(filecmp.cmp('input/gmd_test.txt','input/gmd_test_size0.txt'))


class geneRegionTest(unittest.TestCase):

    def getARecord(self):
        vcfr = VcfReader('input/chr11.subsamples.vcf.gz')
        return vcfr.info_rec

    def test_RL_collist(self):
        collist = '2,4,6'
        rl = RegionList('input/gr_ex_multicolumn.txt', colstr=collist, zeroho=True)
        r = rl.regions[0]
        set_r = Region(100000, 100100, '11')
        self.assertEqual(r, set_r)

    def test_RL_sort(self):
        rl = RegionList('input/unsorted_regions.txt')
        sl = RegionList('input/sorted_regions.txt')
        rl.regions.sort()
        self.assertTrue(rl.regions == sl.regions)

    def test_RL_charsort(self):
        rl = RegionList('input/unsorted_regions.txt',sortmethod="string")
        sl = RegionList('input/charsorted_regions.txt')
        rl.regions.sort()
        self.assertTrue(rl.regions == sl.regions)


    def test_RL_overlap(self):
        rl = RegionList('input/overlap_regions.txt')
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


    def test_generateSequence_insert(self):
        vcf_to_seq(['input/chr11.subsamples.vcf.gz',
             'input/human_g1k_chr11.fasta',
             'input/insert_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('input/chr11.subsamples.fasta',
                         'input/chr11.insert.example.fasta'), True)

    def test_generateSequence_del(self):
        vcf_to_seq(['input/chr11.subsamples.vcf.gz',
             'input/human_g1k_chr11.fasta',
             'input/del_region.txt', '--indels'])
        self.assertEqual(filecmp.cmp('input/chr11.subsamples.fasta',
                         'input/chr11.del.example.fasta'), True)

    def test_generateSequence_multi(self):
        vcf_to_seq(['input/chr11.subsamples.vcf.gz',
             'input/human_g1k_chr11.fasta',
             'input/multi_inregion.txt'])
        self.assertEqual(filecmp.cmp('input/chr11.subsamples.fasta',
                         'input/chr11.multi.example.fasta'), True)

    def test_subsample(self):
        vcf_to_seq(['input/chr11.vcf.gz',
             'input/human_g1k_chr11.fasta',
             'input/snp_region.txt',
             '--subsamp-list', 'input/subsample_list.txt'])
        self.assertEqual(filecmp.cmp('input/chr11.fasta',
                         'input/chr11.snp.example.fasta'), True)

    def test_generateSequence_unzipped(self):
        vcf_to_seq(['input/chr11.unzipped.vcf',
             'input/human_g1k_chr11.fasta',
             'input/snp_region.txt'])
        self.assertEqual(filecmp.cmp('input/chr11.unzipped.fasta',
                         'input/chr11.snp.example.fasta'), True)

    def test_generateSequence_uz_adj(self):
        vcf_to_seq(['input/chr11.unzipped.vcf',
             'input/human_g1k_chr11.fasta',
            'input/adj_regions.txt'])
        self.assertEqual(filecmp.cmp('input/chr11.unzipped.fasta',
                         'input/chr11.adjacent.fasta'), True)


    def test_generateSequence_mismatchref(self):
        args = ['input/chr11.subsamples.vcf.gz',
             'input/human_g1k_messy.fasta',
             'input/snp_region.txt']

        self.assertRaises(Exception, vcf_to_seq, args)

    def tearDown(self):
        tryRemove('input/chr11.subsamples.fasta')
        tryRemove('input/chr11.fasta')
        tryRemove('input/chr11.unzipped.fasta')

class checkCompressionTest(unittest.TestCase):

    def test_checkFormat_nozip(self):
        comp = checkFormat('input/chr11.unzipped.vcf')
        self.assertTrue(comp == 'vcf')

    def test_checkFormat_bgzip(self):
        comp = checkFormat('input/chr11.subsamples.vcf.gz')
        self.assertTrue(comp == 'bgzip')

    def test_checkFormat_gzip(self):
        comp = checkFormat('input/chr11.gzipped.vcf.gz')
        self.assertTrue(comp=='gzip')

    def test_checkFormat_bcf(self):
        comp = checkFormat('input/chr11.test.bcf')
        self.assertTrue(comp=='bcf')

    def test_checkFormat_other(self):
        comp = checkFormat('input/plain.txt')
        self.assertTrue(comp=="other")

    def test_checkFormat_otherzip(self):
        comp = checkFormat('input/plain.txt.gz')
        self.assertTrue(comp=="other")

class reduceTest(unittest.TestCase):

    def test_vcf_region_write_simple(self):
        vcf_region_write(['input/chr11.subsamples.vcf.gz',
                         '--bed','input/snp_region.txt',
                         '--out','input/chr11.test.vcf.gz'])
        self.assertTrue(filecmp.cmp('input/chr11.snp.sr.vcf.gz',
                        'input/chr11.test.vcf.gz'))

    def test_vcf_region_write_unzipped(self):
        vcf_region_write(['input/chr11.unzipped.vcf',
                          '--bed','input/snp_region.txt',
                          '--out','input/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('input/chr11.snp.sr.vcf',
                        'input/chr11.test.vcf'))

    def test_vcf_region_write_first(self):
        vcf_region_write(['input/chr11.unzipped.vcf',
                          '--bed','input/first_region.txt',
                          '--out','input/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('input/chr11.first.vcf',
                        'input/chr11.test.vcf'))

    def test_vcf_region_write_last(self):
        vcf_region_write(['input/chr11.unzipped.vcf',
                          '--bed','input/last_region.txt',
                          '--out','input/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('input/chr11.last.vcf',
                        'input/chr11.test.vcf'))

    def test_vcf_region_write_multi(self):
        vcf_region_write(['input/chr11.unzipped.vcf',
                          '--bed','input/multi_big_region.txt',
                          '--out','input/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('input/chr11.multi.big.vcf',
                        'input/chr11.test.vcf'))

    def test_vcf_region_write_adj(self):
        vcf_region_write(['input/chr11.unzipped.vcf',
                          '--bed','input/adj_big_region.txt',
                          '--out','input/chr11.test.vcf'])
        self.assertTrue(filecmp.cmp('input/chr11.adj.big.vcf',
                        'input/chr11.test.vcf'))

    def tearDown(self):
        pass
        #tryRemove('input/chr11.test.vcf.gz')
        #tryRemove('input/chr11.test.vcf')



class configTest(unittest.TestCase):

    def test_vcf_to_seq_config(self):
        vcf_to_seq(['input/chr11.subsamples.vcf.gz',
                    'input/human_g1k_chr11.fasta',
                    'input/insert_region.txt',
                    '--conf','input/config.config'])
        self.assertTrue(filecmp.cmp('input/chr11.insert.example.fasta',
                        'input/chr11.subsamples.fasta'))

    def test_vcf_to_seq_config_required(self):
        vcf_to_seq(['--conf', 'input/withreq.conf'])
        self.assertTrue(filecmp.cmp('input/chr11.snp.example.fasta',
                        'input/chr11.subsamples.fasta'))

class fourgameteTest(unittest.TestCase):

    def test_hk(self):
        region_list = [[[192385, 196945], [198986, 202254], [202547, 202813],
         [206421, 207401], [207400, 215365], [217814, 218641]]]
        test_list = sample_fourgametetest_intervals(["--vcfs","input/chr11subsamples4gtest.vcf.gz","--hk","--retl"])
        self.assertEqual(region_list, test_list)

    def test_comp(self):
        region_list = [[[192385, 196943], [192386, 202252], [198987, 202546],
         [201585, 202811], [202548, 207399], [206422, 215363], [207401, 218140],
          [207463, 218639], [217815, 221584]]]
        test_list = sample_fourgametetest_intervals(["--vcfs","input/chr11subsamples4gtest.vcf.gz","--4gcompat","--retl",'--ovlpi','--ovlps'])
        self.assertEqual(region_list,test_list)

    def test_vcf(self):
        sample_fourgametetest_intervals(["--vcfs","input/chr11.subsamples.vcf.gz",'--4gcompat','--reti','--out','input/chr11.4g.vcf','--numinf','2'])
        self.assertTrue(filecmp.cmp('input/chr11.4gtest1.vcf','input/chr11.4g.vcf'))
        tryRemove('input/chr11.4g.vcf')

    def test_consecutive(self):
        sample_fourgametetest_intervals(["--vcfs","input/chr11.4gsmall.vcf",'--4gcompat','--reti','--out','input/chr11.4g.vcf','--numinf','2','--ovlpi','--ovlps'])
        self.assertTrue(filecmp.cmp('input/chr11.4gtest2.vcf','input/chr11.4g.vcf'))
        tryRemove('input/chr11.4g.vcf')


class imaTest(unittest.TestCase):
    def test_ima(self):
        vcf_to_ima(['--vcf','input/chr11.subsamples.vcf.gz',
                    '--reference-fasta','input/human_g1k_chr11.fasta',
                    '--bed','input/snp_region.txt',
                    '--model-file','input/testmodel.model',
                    '--zero-ho',
                    '--out','input/chr11.subsamples.ima.u'])
        self.assertTrue(filecmp.cmp('input/chr11.ima.u',
                        'input/chr11.subsamples.ima.u'))

class phaseTest(unittest.TestCase):
    def test_shapeit(self):
        vcf_phase.run(['--vcf','input/chr11.adj.big.vcf',
                   '--out','input/chr11.test.vcf',
                   '--out-format','vcf',
                   '--phase-algorithm','shapeit',
                   '--random-seed','123','--overwrite'])
        self.assertTrue(filecmp.cmp('input/chr11.test.vcf',
                        'input/chr11.shapeit.test.vcf'))

    def test_beagle(self):
        vcf_phase.run(['--vcf','input/chr11.adj.big.vcf',
                   '--out','input/chr11.test.vcf',
                   '--out-format','vcf',
                   '--phase-algorithm','beagle',
                   '--random-seed','123','--overwrite'])
        self.assertTrue(filecmp.cmp('input/chr11.test.vcf',
                        'input/chr11.beagle.test.vcf'))

    def tearDown(self):
        tryRemove('input/chr11.test.vcf')



if __name__ == "__main__":
    initLogger(filelevel="ERROR",streamlevel="ERROR")
    unittest.main(verbosity=2)
