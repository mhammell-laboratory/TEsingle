import copy
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from TEtools.Counters import UMI_CATEGORIES, basic_filtering, count_cells
from TEtools.GeneIndex import GFF_Reader, GeneFeatures
from TEtools.TEindex import TEfeatures


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_gene_ids_are_unquoted_with_gzip_gtf(self):
        path = self.root / 'genes.gtf.gz'
        with gzip.open(path, 'wt') as stream:
            stream.write('1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "GENE1";\n\n')
        self.assertEqual([row[0] for row in GFF_Reader(str(path), 'gene_id')], ['GENE1'])

    def test_gene_feature_order_and_transcript_intervals_are_preserved(self):
        path = self.root / 'genes.gtf'
        path.write_text(''.join('1\ttest\t%s\t%s\t%s\t.\t+\t.\tgene_id "%s";\n' % row for row in [
            ('transcript', 1, 100, 'B'), ('transcript', 200, 300, 'A'),
            ('transcript', 400, 500, 'B'), ('exon', 600, 700, 'IGNORED')]))
        index = GeneFeatures(str(path), 'forward', 'transcript', 'gene_id')
        self.assertEqual(index.getFeatures(), ['B', 'A'])
        self.assertEqual(index.Gene_annotation([['1', 450, 460, '+']]), [0])

    def test_TE_attributes_allow_standard_gtf_whitespace(self):
        path = self.root / 'tes.gtf'
        path.write_text('1\ttest\texon\t1\t100\t.\t+\t.\t'
                        'gene_id "L1";  transcript_id "L1-1"; family_id "L1"; class_id "LINE" ;\n\n')
        index = TEfeatures()
        index.build(str(path))
        self.assertEqual(index.getElements(), ['L1-1:L1:L1:LINE'])
        self.assertEqual(index.TE_annotation([['1', 1, 50, '+']]), [0])


class CountingTests(unittest.TestCase):
    def test_sequential_counts_keep_cells_without_umis(self):
        cells = {'CELL2': {'AAAAAAAAAAAA': {'read': [([0], [])]}}, 'CELL1': {}}
        filtered, before, after = basic_filtering(copy.deepcopy(cells), 0)
        self.assertEqual((before, after), (2, 2))
        matrix, barcodes, _ = count_cells(filtered, sparse.lil_matrix((1, 1)), sparse.lil_matrix((1, 1)), 1)
        self.assertEqual(barcodes, ['CELL1', 'CELL2'])
        np.testing.assert_array_equal(matrix.toarray(), [[0, 0], [1, 0]])

    def test_umi_categories_separate_TE_only_and_unannotated(self):
        # UMIs differ at six or more positions, so none are merged.
        # Each read maps to a list of (gene indexes, TE indexes) alignments.
        cells = {'CELL_A': {'AAAAAAAAAAAA': {'te_multi': [([], [0]), ([], [1])]},
                            'CCCCCCCCCCCC': {'mixed_multi': [([0], []), ([], [1])]},
                            'GGGGGGGGGGGG': {'gene_multi': [([0], []), ([1], [])]},
                            'TTTTTTTTTTTT': {'te_unique': [([], [1])]},
                            'ACACACACACAC': {'both_unique': [([0], [0])]},
                            'GTGTGTGTGTGT': {'gene_unique': [([1], [])]},
                            'CACACACACACA': {'intergenic_unique': [([], [])],
                                             'intergenic_multi': [([], []), ([], [])]}},
                 'CELL_B': {'AAAAAAAAAAAA': {'te_multi': [([], [0]), ([], [1])]},
                            'CCCCCCCCCCCC': {'intergenic_unique': [([], [])]}}}
        matrix, barcodes, summary = count_cells(cells, sparse.lil_matrix((1, 2)), sparse.lil_matrix((1, 2)), 1)
        self.assertEqual(dict(zip(UMI_CATEGORIES, summary.tolist())),
                         {'unmoored_TE': 2, 'unmoored_gene_TE': 1, 'unmoored_gene': 1, 'anchored_TE': 1,
                          'anchored_gene_TE': 1, 'anchored_gene': 1, 'unannotated': 2})
        counts = dict(zip(barcodes, matrix.toarray()))
        self.assertAlmostEqual(counts['CELL_A'].sum(), 6)
        self.assertGreaterEqual(counts['CELL_A'][2:].sum(), 2)
        np.testing.assert_array_equal(counts['CELL_B'][:2], [0, 0])
        self.assertAlmostEqual(counts['CELL_B'][2:].sum(), 1)


if __name__ == '__main__':
    unittest.main()
