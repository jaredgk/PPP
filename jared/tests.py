import unittest
from vcf_ref_to_seq import getMaxAlleleLength



class FirstTestOfGmal(unittest.TestCase):
    def runTest(self):
        l = ('AAA','A','AAAAAAA')
        self.assertEqual(getMaxAlleleLength(l),7)

if __name__ == "__main__":
    unittest.main()
