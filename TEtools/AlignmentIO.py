import pysam


def annotation_direction(stranded, read_is_reverse):
    if stranded == 'no':
        return '.'
    else:
        direction = [1, -1][stranded == 'reverse']
        algn_dir = [1, -1][read_is_reverse]
        return ['.', '+', '-'][algn_dir * direction]


def read_in_alignment(filename, geneIdx, teIdx, stranded):
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

    return transc_dict
