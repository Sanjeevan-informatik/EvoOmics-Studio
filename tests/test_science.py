import math
from pathlib import Path
import tempfile
import unittest

from evoomics.io.sequences import read_fasta, read_vcf
from evoomics.models import Alignment, ValidationError
from evoomics.phylogeny.distance import pair_distance, distance_matrix
from evoomics.phylogeny.tree import upgma
from evoomics.qc import filter_alignment

ROOT = Path(__file__).resolve().parents[1]


class AlignmentTests(unittest.TestCase):
    def test_invalid_alignments(self):
        cases = [(('a', 'a'), ('AA', 'AA')), (('a', 'b'), ('A', 'AA')),
                 (('a', 'b'), ('', '')), (('a', 'b'), ('AU', 'AA')),
                 (('a b', 'c'), ('AA', 'AA')), (('a',), ('AA',))]
        for names, seqs in cases:
            with self.subTest(names=names, seqs=seqs), self.assertRaises(ValidationError):
                Alignment(names, seqs)

    def test_filter(self):
        result, qc, keep = filter_alignment(Alignment(('a', 'b'), ('ACN', 'ATG')), 0)
        self.assertEqual(result.sequences, ('AC', 'AT'))
        self.assertEqual(keep, [0, 1])
        self.assertEqual(qc['removed_sites'], 1)

    def test_filter_empty_and_invalid(self):
        for threshold in (0, -0.1, 1.1, float('nan')):
            with self.subTest(threshold=threshold), self.assertRaises(ValidationError):
                filter_alignment(Alignment(('a', 'b'), ('N', 'N')), threshold)


class DistanceTests(unittest.TestCase):
    def test_known_p(self):
        self.assertEqual(pair_distance('ACGTN-R', 'ATGTAAA'), (0.25, 4))

    def test_jc69_hand_calculation(self):
        d, n = pair_distance('AAAA', 'AAAG', 'jc69')
        self.assertAlmostEqual(d, -0.75 * math.log(2/3))
        self.assertEqual(n, 4)

    def test_k2p_hand_calculation(self):
        d, n = pair_distance('AAAAAAAAAA', 'GCAAAAAAAA', 'k2p')
        self.assertAlmostEqual(d, -0.5 * math.log(0.7) - 0.25 * math.log(0.8))
        self.assertEqual(n, 10)

    def test_saturation(self):
        for model in ('jc69', 'k2p'):
            with self.subTest(model=model), self.assertRaises(ValidationError):
                pair_distance('AAAA', 'GGGG', model)

    def test_insufficient_overlap(self):
        for a, b, minimum in [('NN', 'AA', 1), ('AN', 'AA', 2)]:
            with self.subTest(a=a), self.assertRaises(ValidationError):
                pair_distance(a, b, min_sites=minimum)

    def test_matrix(self):
        matrix, counts = distance_matrix(Alignment(('a', 'b'), ('AN', 'AT')))
        self.assertEqual(matrix, [[0, 0], [0, 0]])
        self.assertEqual(counts, [[1, 1], [1, 2]])


class TreeTests(unittest.TestCase):
    def test_expected_tree(self):
        self.assertEqual(upgma(['a','b','c'], [[0,2,4],[2,0,4],[4,4,0]]),
                         "('c':2,('a':1,'b':1):1);")

    def test_size_weighted_merge(self):
        tree = upgma(['a','b','c','d'], [[0,2,4,8],[2,0,4,8],[4,4,0,10],[8,8,10,0]])
        self.assertEqual(tree, "('d':4.333333333,('c':2,('a':1,'b':1):1):2.333333333);")

    def test_invalid_matrices(self):
        for matrix in ([[0,1],[2,0]], [[0,-1],[-1,0]], [[0,float('inf')],[float('inf'),0]], [[1,0],[0,0]]):
            with self.subTest(matrix=matrix), self.assertRaises(ValidationError):
                upgma(['a','b'], matrix)


class ReaderTests(unittest.TestCase):
    def load_text(self, text, reader):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'input'
            p.write_text(text)
            return reader(p)

    def test_fasta(self):
        alignment, _ = self.load_text('>a description\nac\ngt\n>b\nACGT\n', read_fasta)
        self.assertEqual(alignment.sequences, ('ACGT', 'ACGT'))

    def test_bad_fasta(self):
        for text in ('AA\n>a\nAA', '>\nAA', '>a\nAA\n>a\nAA'):
            with self.subTest(text=text), self.assertRaises(ValidationError):
                self.load_text(text, read_fasta)

    def test_vcf_mask_and_iupac(self):
        alignment, stats = read_vcf(ROOT / 'examples/tiny.vcf')
        self.assertEqual(alignment.sequences, ('ACGT', 'ATNT', 'GTGT'))
        self.assertEqual(stats['heterozygous_genotypes'], 1)
        alignment, _ = read_vcf(ROOT / 'examples/tiny.vcf', 'iupac')
        self.assertEqual(alignment.sequences[1], 'ATRT')

    def test_vcf_missing_phased_haploid_and_gt_position(self):
        text = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n'
        text += 'chr1\t1\t.\tA\tG\t.\tPASS\t.\tDP:GT\t9:1|1\t8:0\t7:0/.\n'
        a, stats = self.load_text(text, read_vcf)
        self.assertEqual(a.sequences, ('G','A','N'))
        self.assertEqual(stats['missing_genotypes'], 1)

    def test_vcf_invalid_and_duplicate(self):
        text = (ROOT / 'examples/tiny.vcf').read_text()
        for modified in (text.replace('0/1', '0/2'), text.replace('0/1', '0/0/1'),
                         text + text.splitlines()[-1] + '\n', text.replace('FORMAT', 'BAD')):
            with self.subTest(text=modified), self.assertRaises(ValidationError):
                self.load_text(modified, read_vcf)

    def test_vcf_skips_filtered_indel_multiallelic(self):
        text = (ROOT / 'examples/tiny.vcf').read_text()
        row = text.splitlines()[-1]
        text += row.replace('40', '50').replace('PASS', '.') + '\n'
        text += row.replace('40', '60').replace('\tC\t', '\tCA\t') + '\n'
        text += row.replace('40', '70').replace('\tC\t', '\tC,G\t') + '\n'
        a, stats = self.load_text(text, read_vcf)
        self.assertEqual(a.sites, 4)
        self.assertEqual(stats['skipped_filter'], 1)
        self.assertEqual(stats['skipped_non_biallelic_snp'], 2)
