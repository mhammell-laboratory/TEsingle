import scipy.sparse as sparse
import numpy as np
import multiprocessing as mp
import time

from TEtools.CellGrapher import *
from TEtools.EM import *


def normalizelil(sparse_vect):
    tmp_sum = np.sum(sparse_vect.data[0])  # has to take 0th index for lil structure
    if tmp_sum > 0:
        norm = 1.0 / tmp_sum
    else:
        norm = 0.0
    return sparse_vect.multiply(norm)


def rearranger(subg_list):
    transc_dict = {'gene_only': [], 'tes': [], 'anch_both': []}
    cat_sum = np.zeros(6)

    for subg_assgn in subg_list:
        scnt = subg_assgn['cnt']
        cat_sum += subg_assgn['cat']
        if np.sum(subg_assgn['cat']*[0,0,1,0,0,1]) > 0:  # gene-only, unmoored or anchored
            transc_dict['gene_only'].append(scnt)
        elif np.sum(subg_assgn['cat']*[0,0,0,0,1,0]) > 0:  # anchored, gene+TE
            transc_dict['anch_both'].append(scnt)
        else:
            transc_dict['tes'].append(scnt)

    return transc_dict['gene_only'], transc_dict['tes'], transc_dict['anch_both'], cat_sum


def resolver(categories, both_list):
    denom = categories[3] + categories[5]
    if denom == 0:
        te_frac = 0
    else:
        te_frac = categories[3]/denom
    gene_frac = 1-te_frac
    sum_list = []
    for umi in both_list:
        count = umi[0].multiply(gene_frac)+umi[1].multiply(te_frac)
        count = count.tolil()
        count = normalizelil(count)
        sum_list.append(count)

    return sum_list


def summer(cnt_list, blank):
    tot_cnt = blank.copy()
    for vect in cnt_list:
        tot_cnt += vect

    tot_cnt = sparse.csr_matrix(tot_cnt)
    return tot_cnt


def count_reads(gene_blank, te_blank, cbc_dict, nodelist):
    category = np.zeros(6)
    anchored = []
    equiv_classes = {}

    for umi in nodelist:
        for read in cbc_dict[umi]:
            equiv_classes[read] = {'genes': [], 'tes': []}

            if len(cbc_dict[umi][read]) == 1:
                anchored.append(read)

            for alignment in cbc_dict[umi][read]:
                tmp_genes = alignment[0]
                tmp_tes = alignment[1]

                if len(tmp_genes) > 0:
                    equiv_classes[read]['genes'].append(tmp_genes)  # this will be a list of lists of genes

                if len(tmp_tes) > 0:
                    equiv_classes[read]['tes'].append(tmp_tes)
        del cbc_dict[umi]

    for read in anchored[:]:  # remove any reads that uniquely mapped but didn't overlap any annotations
        if len(equiv_classes[read]['genes']) == 0 and len(equiv_classes[read]['tes']) == 0:
            anchored.remove(read)

    if len(anchored) == 0:  # unmoored umi
        unmo_gene_cnts = gene_blank.copy()
        unmo_te_cnts = te_blank.copy()
        tes = False
        for read in equiv_classes:
            num_te_algn = len(equiv_classes[read]['tes'])
            num_ge_algn = len(equiv_classes[read]['genes'])

            if num_te_algn > 0:
                tes = True
                for algn in equiv_classes[read]['tes']:
                    num_te_annot = len(algn)
                    for annot in algn:
                        unmo_te_cnts[0,annot] += 1 / (num_te_algn * num_te_annot)

            if num_ge_algn > 0:
                for algn in equiv_classes[read]['genes']:
                    num_ge_annot = len(algn)
                    for annot in algn:
                        unmo_gene_cnts[0,annot] += 1 / (num_ge_algn * num_ge_annot)

        if tes:
            unmo_te_cnts = normalizelil(unmo_te_cnts)
            if np.sum(unmo_gene_cnts.data) == 0:
                category[0] = 1
                unmo_cnts = sparse.hstack((gene_blank, unmo_te_cnts))
            else:
                category[1] = 1
                unmo_gene_cnts = normalizelil(unmo_gene_cnts)
                unmo_cnts = sparse.hstack((unmo_gene_cnts.multiply(0.5), unmo_te_cnts.multiply(0.5)))

            return {'cnt': unmo_cnts, 'cat': category}

        else:
            category[2] = 1
            unmo_gene_cnts = normalizelil(unmo_gene_cnts)
            return {'cnt': sparse.hstack((unmo_gene_cnts, te_blank)), 'cat': category}

    else:  # anchored umi - consider making intersection possible instead of just union?
        anchored_annotations = [set(), set()]
        extras = False
        for read in anchored:  # anchored read has only one alignment, but could have both gene and te annotations
            if len(equiv_classes[read]['genes']) > 0:
                anchored_annotations[0] = anchored_annotations[0].union(set(equiv_classes[read]['genes'][0]))
            if len(equiv_classes[read]['tes']) > 0:
                anchored_annotations[1] = anchored_annotations[1].union(set(equiv_classes[read]['tes'][0]))
        if len(anchored_annotations[0]) > 0:
            if len(anchored_annotations[1]) > 0:
                category[4] = 1
                extras = True
                anch_cnts_extra = te_blank.copy()
                anchored_annotations_extra = anchored_annotations[1]
            else:
                category[5] = 1
            anchored_annotations = anchored_annotations[0]
            anch_cnts = gene_blank.copy()

        else:  # anchored reads are required to have annotations, so if no genes then grab tes
            anchored_annotations = anchored_annotations[1]
            anch_cnts = te_blank.copy()
            category[3] = 1

        anchored_count_dict = dict.fromkeys(anchored_annotations, 0)
        for read in equiv_classes:
            for alignment in equiv_classes[read][['genes', 'tes'][int(category[3])]]:
                for val in anchored_annotations.intersection(set(alignment)):
                    anchored_count_dict[val] += 1

        for val in anchored_count_dict:
            anch_cnts[0,val] = anchored_count_dict[val]

        anch_cnts = normalizelil(anch_cnts)
        if category[3] != 1:
            anch_cnts = sparse.hstack((anch_cnts, te_blank))
        else:
            anch_cnts = sparse.hstack((gene_blank, anch_cnts))

        if extras:
            anchored_cnt_dict_extra = dict.fromkeys(anchored_annotations_extra, 0)
            for read in equiv_classes:
                for alignment in equiv_classes[read]['tes']:
                    for val in anchored_annotations_extra.intersection(set(alignment)):
                        anchored_cnt_dict_extra[val] += 1

            for val in anchored_cnt_dict_extra:
                anch_cnts_extra[0,val] = anchored_cnt_dict_extra[val]

            anch_cnts_extra = normalizelil(anch_cnts_extra)
            anch_cnts_extra = sparse.hstack((gene_blank, anch_cnts_extra))

            return {'cnt': [anch_cnts, anch_cnts_extra], 'cat': category}

        return {'cnt': anch_cnts, 'cat': category}


def count_umis(cell_dict, gene_blank, te_blank):
    bc_graph = cell_graph()
    bc_graph.build(cell_dict)

    subg_counts = []

    for subgraph in bc_graph.subgraphs:
        subg_counts.append(count_reads(gene_blank, te_blank, cell_dict, list(subgraph.nodes)))

    del bc_graph

    gene_vects, te_vects, anchboth_vects, cat_sum = rearranger(subg_counts)
    anchboth_vects = resolver(cat_sum, anchboth_vects)
    gene_cnts = summer(gene_vects, sparse.hstack((gene_blank, te_blank)))

    EM_vects = te_vects + anchboth_vects

    tot_cnt = EM(EM_vects, sparse.hstack((gene_blank, te_blank))) + gene_cnts
    return tot_cnt, cat_sum


def count_cells(full_dict, gene_blank, te_blank, num_processes):
    cbclist = list(full_dict.keys())

    mp.set_start_method('spawn')

    pool = mp.Pool(num_processes)

    tbl_list = [(cbc, pool.apply_async(count_umis, args=(full_dict[cbc], gene_blank, te_blank))) for cbc in cbclist]
    del full_dict
    full_results = [(cbc, n.get()) for (cbc, n) in tbl_list]

    pool.close()
    pool.join()

    print('count_umis finishing at:  ', time.ctime(), flush=True)
    checker = False

    new_cbclist = []
    full_sum = np.zeros(6)
    tblstack = False

    full_results.sort(key=lambda x: x[0])

    for (cbc,placeholder) in full_results:
        cell_cnt = placeholder[0]
        cat_sum = placeholder[1]
        new_cbclist.append(cbc)
        full_sum += cat_sum

        if not checker:
            tblstack = cell_cnt
            checker = True
        else:
            tblstack = sparse.vstack([tblstack, cell_cnt])

    print('count_cells finishing at:  ', time.ctime(), flush=True)
    return tblstack, new_cbclist, full_sum


def basic_filtering(bcumi_dict, cutoff):
    init_cbc = len(bcumi_dict)
    for cbc in bcumi_dict.copy():
        if len(bcumi_dict[cbc]) < cutoff:
            del bcumi_dict[cbc]
    post_cbc = len(bcumi_dict)
    return bcumi_dict, init_cbc, post_cbc