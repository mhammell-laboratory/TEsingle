import gzip
import re
import shlex

import pysam


def normalize_barcode(barcode):
    """Remove the 10x GEM suffix without changing the nucleotide barcode."""
    return re.sub(r'-\d+$', '', barcode.strip())


def read_cell_barcodes(filename):
    opener = gzip.open if filename.endswith('.gz') else open
    with opener(filename, 'rt') as stream:
        barcodes = {normalize_barcode(line.split()[0]) for line in stream if line.strip()}
    if not barcodes:
        raise ValueError('Cell barcode list is empty: ' + filename)
    return barcodes


def validate_alignment(filename, allow_truncated=False):
    """Catch the STAR setting that makes a multimapper look uniquely aligned."""
    with pysam.AlignmentFile(filename, check_sq=False) as bam:
        for program in bam.header.to_dict().get('PG', []):
            if program.get('PN') != 'STAR':
                continue
            tokens = shlex.split(program.get('CL', ''))
            def option(name, default):
                return int(tokens[tokens.index(name) + 1]) if name in tokens else default
            reported = option('--outSAMmultNmax', -1)
            permitted = option('--outFilterMultimapNmax', 10)
            if 0 < reported < permitted and not allow_truncated:
                raise ValueError(
                    'STAR BAM has truncated multimapping alignments '
                    '(--outSAMmultNmax %s, --outFilterMultimapNmax %s). '
                    'Realign FASTQs with --outSAMmultNmax -1. '
                    'The missing alignments cannot be recovered from this BAM. '
                    '--allow-truncated-bam is for diagnostic comparisons only.' % (reported, permitted))


def annotation_direction(stranded, read_is_reverse):
    if stranded == 'no':
        return '.'
    else:
        direction = [1, -1][stranded == 'reverse']
        algn_dir = [1, -1][read_is_reverse]
        return ['.', '+', '-'][algn_dir * direction]


def read_in_alignment(filename, geneIdx, teIdx, stranded, cell_barcodes=None):
    # given .sam or .bam file, returns dictionary of {barcodes:UMIs} and {UMIs:{reads:{alignments}}}
    # alignments are in format [chromosome, [0-indexed blocks in tuples], direction]

    psave = pysam.set_verbosity(0)  # disable htslib verbosity to avoid "no index file" warning
    samfile = pysam.AlignmentFile(filename, check_sq=False)
    pysam.set_verbosity(psave)  # re-enable htslib verbosity

    transc_dict = {}
    err_reads = 0

    try:
        while 1:
            aligned_read = next(samfile)
            if aligned_read.is_unmapped or aligned_read.is_duplicate or aligned_read.is_qcfail:
                continue

            rbc = False
            rumi = False

            try:
                rbc = aligned_read.get_tag('CB')
                if rbc == '-':
                    rbc = False
            except:
                err_reads += 1

            if not rbc:
                continue

            rbc = normalize_barcode(rbc)
            if cell_barcodes is not None and rbc not in cell_barcodes:
                continue

            try:
                rumi = aligned_read.get_tag('UR')
                if rumi == '-':
                    rumi = False

            except:
                err_reads += 1

            if not (rbc and rumi):
                continue

            if rbc not in transc_dict:
                transc_dict[rbc] = {}

            if rumi not in transc_dict[rbc]:
                transc_dict[rbc][rumi] = {}

            cur_read_name = aligned_read.query_name

            if cur_read_name not in transc_dict[rbc][rumi]:
                transc_dict[rbc][rumi][cur_read_name] = []

            algn = aligned_read.get_blocks()  # 0-indexed
            chrom = aligned_read.reference_name
            direction = annotation_direction(stranded, aligned_read.is_reverse)

            itv_list = []

            for block in algn:
                itv_list.append([chrom, block[0]+1, block[1], direction])

            tmp_genes = geneIdx.Gene_annotation(itv_list)
            tmp_tes = teIdx.TE_annotation(itv_list)
            transc_dict[rbc][rumi][cur_read_name].append(tuple([tmp_genes, tmp_tes]))

    except StopIteration:
        pass

    finally:
        samfile.close()

    return transc_dict
