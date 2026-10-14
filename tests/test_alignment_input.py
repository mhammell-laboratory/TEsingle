import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from TEtools.AlignmentIO import read_cell_barcodes, read_in_alignment, validate_alignment


class FakeIndex:
    def Gene_annotation(self, intervals):
        return [0]

    def TE_annotation(self, intervals):
        return []


class AlignmentInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def make_bam(self, limit, alignments=1):
        path = self.root / 'test.bam'
        header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': '1', 'LN': 10000}],
                  'PG': [{'ID': 'STAR', 'PN': 'STAR',
                          'CL': 'STAR --outSAMmultNmax %s --outFilterMultimapNmax 100' % limit}]}
        with pysam.AlignmentFile(str(path), 'wb', header=header) as stream:
            for i in range(alignments):
                read = pysam.AlignedSegment()
                read.query_name = 'multimapper'
                read.query_sequence = 'A' * 50
                read.reference_id = 0
                read.reference_start = 100 + i * 100
                read.flag = 0 if i == 0 else 256
                read.cigarstring = '50M'
                read.set_tags([('CB', 'ACGTACGTACGTACGT-1'), ('UR', 'AAAAAAAAAAAA'), ('NH', 2)])
                stream.write(read)
        return path

    def test_reject_truncated_multimapping_bam(self):
        with self.assertRaisesRegex(ValueError, 'missing alignments cannot be recovered'):
            validate_alignment(str(self.make_bam(1)))

    def test_accept_all_alignments_and_diagnostic_override(self):
        validate_alignment(str(self.make_bam(-1, 2)))
        validate_alignment(str(self.make_bam(1)), allow_truncated=True)

    def test_called_cell_list_accepts_gzip_and_10x_suffix(self):
        path = self.root / 'cells.tsv.gz'
        with gzip.open(path, 'wt') as stream:
            stream.write('ACGTACGTACGTACGT-1\nACGTACGTACGTACGT\n')
        self.assertEqual(read_cell_barcodes(str(path)), {'ACGTACGTACGTACGT'})

    def test_keep_all_multimapping_records_and_filter_cells(self):
        path = self.make_bam(-1, 2)
        result = read_in_alignment(str(path), FakeIndex(), FakeIndex(), 'forward', {'ACGTACGTACGTACGT'})
        self.assertEqual(len(result['ACGTACGTACGTACGT']['AAAAAAAAAAAA']['multimapper']), 2)
        self.assertEqual(read_in_alignment(str(path), FakeIndex(), FakeIndex(), 'forward', {'OTHER'}), {})


if __name__ == '__main__':
    unittest.main()
