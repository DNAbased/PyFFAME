"""Unittests for PyFFAME Framework"""

import unittest
from handlers import BarcodeHandler, EnzymeHandler, GenomicHandler, DataHelper, MongoHandler


class TestHelperFunctions(unittest.TestCase):  
    
    def test_revcomp(self):
        dh = DataHelper()
        start_dna_seqence = "ATTACGATTACG"
        expected_revcomp_sequence = "CGTAATCGTAAT"
        processed_dna_sequence = dh.revcomp(start_dna_seqence)
        self.assertEqual(processed_dna_sequence, expected_revcomp_sequence)

    def test_rev(self):
        dh = DataHelper()
        start_dna_seqence = "ATTACGATTACG"
        expected_rev_sequence = "GCATTAGCATTA"
        processed_dna_sequence = dh.rev(start_dna_seqence)
        self.assertEqual(processed_dna_sequence, expected_rev_sequence)

    def test_comp(self):
        dh = DataHelper()
        start_dna_seqence = "ATTACGATTACG"
        expected_comp_sequence = "TAATGCTAATGC"
        processed_dna_sequence = dh.comp(start_dna_seqence)
        self.assertEqual(processed_dna_sequence, expected_comp_sequence)


class TestBarcodeHandlerFunctions(unittest.TestCase): 
    def test_first_thing(self):
        pass


if __name__ == '__main__':
    unittest.main()
