TEsingle
=============

Version: 0.1

*NOTE* TEsingle relies on specially curated and indexed GTF files, which are not
packaged with this software due to their size. Please go to
`our website <http://hammelllab.labsites.cshl.edu/software/tesingle>`_
for instructions to download the prebuilt curated indices.

TEsingle takes single-cell RNA-seq data and annotates transcripts to both
genes & transposable elements, producing a count table of all UMI counts for all cell barcodes.


`Github Page <https://github.com/mhammell-laboratory/TEsingle>`_

`Molly Gale Hammell Lab <https://www.mghlab.org/software>`_

Created by Talitha Forcier, Cole Wunderlich, Oliver Tam & Molly Gale Hammell, March 2024

Copyright (C) 2024 Talitha Forcier, Cole Wunderlich, Oliver Tam & Molly Gale Hammell

Contact: mghcompbio@gmail.com

Requirements
------------

Python:     3.2.x or greater

pysam:      0.9.x or greater

networkx

scipy

numpy


Installation
------------

1. Download compressed tarball.
2. Unpack tarball.
3. Navigate into unpacked directory.
4. Run the following::

    $ python setup.py install

If you want to install locally (e.g. /local/home/usr),
run this command instead::

    $ python setup.py install --prefix /local/home/usr

*NOTE* In the above example, you must add
::

    /local/home/usr/bin

to the ``PATH`` variable, and
::

     /local/home/usr/lib/pythonX.Y/site-packages

to the ``PYTHONPATH`` variable, where ``X`` refers to the major python version,
and ``Y`` refers to the minor python version. (e.g. ``python2.7`` if using python version 2.7.x,
and ``python3.6`` if using python version 3.6.x)


TEsingle
=======

Usage
-----

::

    usage: TEsingle -b alignment-file
                   --GTF genic-annot-file
                   --TE TE-annot-file
                   [optional arguments]

    Required arguments:
      -b | --BAM alignment-file    RNAseq alignment file (BAM preferred)
      --GTF genic-annot-file       GTF file for gene annotations
      --TE TE-annot-file           GTF file for transposable element annotations

    Optional arguments:

      *Input/Output options*
      --stranded [option]   Is this a stranded library? (no, forward, or reverse).
                 no      -  Library is unstranded
                 forward -  "Second-strand cDNA library (e.g. 10x Genomics)
                 reverse -  "First-strand" cDNA library (e.g. Illumina TruSeq stranded)
                            DEFAULT: forward.
      --project [name]      Prefix used for output files (e.g. project name)
                            DEFAULT: TEsingle_out

      *Analysis options*
      --cutoff [number]     Minimum number of uncorrected UMIs required to process a barcode
                            DEFAULT: 1000
      --nuc_seq             Flag to set if this is a nuc-seq experiment.
      --locus               Flag for locus-specific run.

      *Other options*
       --threads [number] Number of processors/threads allocated. DEFAULT:10
      -h | --help
         Show help message
      --version
         Show program's version and exit


Example Command Lines
---------------------

::

    TEsingle  --threads 10 --stranded forward -b RNAseq.bam --nuc_seq --GTF refseq_genes.gtf --TE rmsk_TE.gtf --project sample_test

Cluster Usage Recommendations
-----------------------------

In our experience, we recommend around 200-300Gb of memory for analyzing human samples (hg38)
with around 20-30 million mapped reads, when running on a cluster with 10 processors allocated.


Recommendations for TEsingle input files
=============================================

TEsingle can perform transposable element quantification from alignment results (e.g. BAM files)
generated from a variety of programs. Given the variety of experimental systems,
we could not provide an optimal alignment strategy for every approach.
Therefore, we recommend that users identify the optimal parameters for their particular genome
and alignment program in order to get the best results.

When optimizing the alignment parameters, we recommend taking these points into consideration:

*Allowing sufficient number of multi-mappers during alignment*

Most alignment programs provide only 1 alignment per read by default. We recommend reporting multiple
alignments per read. We have found that reporting a maximum of 100 alignments per read provides an
optimal compromise between the size of the alignment file and recovery of multi-mappers in many genome builds.
However, we highly suggest that users optimize this parameter for their particular experiment,
as this could significantly improve the quality of transposable element quantification.

*Specific recommendations when using STAR*

`STAR <https://github.com/alexdobin/STAR>`_ utilizes two parameters for optimal identification of multi-mappers
``--outFilterMultimapNmax`` and ``--winAnchorMultimapNmax``. The author of STAR recommends that ``--winAnchorMultimapNmax``
should be set at twice the value used in ``--outFilterMultimapNmax``, but no less than 50. In our study, we used
100 for ``--outFilterMultimapNmax`` and 200 for ``--winAnchorMultimapNmax``, though we highly suggest users test
multiple values to identify the optimal value for their experiment.

STAR  settings used::

--alignIntronMax 1000000
--alignIntronMin 20
--alignMatesGapMax 1000000
--alignSJDBoverhangMin 1
--alignSJoverhangMin 8
--outFilterMismatchNmax 999
--outFilterMismatchNoverReadLmax 0.04
--outFilterMultimapNmax 100
--winAnchorMultimapNmax 200
--outFilterType BySJout
--outSAMattributes NH HI AS nM CR CY UR UY CB GX GN sS sQ sM
--outSAMheaderHD @HD VN:1.4 SO:unsorted
--outSAMstrandField intronMotif
--outSAMtype BAM SortedByCoordinate
--sjdbScore 1
--soloType CB_UMI_Simple
--soloCellFilter EmptyDrops_CR 5000 0.99 10 45000 90000 500 0.01 20000 0.01 10000
--soloCBmatchWLtype 1MM
--soloUMIdedup 1MM_All
--soloUMIfiltering -
--soloMultiMappers Unique
--soloFeatures GeneFull


Copying & distribution
======================

TEsingle is part of `TEToolkit suite <https://www.mghlab.org/software>`_.

TEsingle is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TEsingle.  If not, see `this website <http://www.gnu.org/licenses/>`_.


