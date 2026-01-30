TEsingle
=============

Version: 1.0

*NOTE* TEsingle relies on specially curated GTF files, which are not
packaged with this software due to their size. Please go to
`our website <https://www.mghlab.org/software/tesingle>`_
for instructions to download the curated annotation files, or
they can be found on `Dropbox <https://www.dropbox.com/scl/fo/jdpgn6fl8ngd3th3zebap/ACdZkShDC1au-OckIipI5kM/TEtranscripts/TE_GTF?rlkey=41oz6ppggy82uha5i3yo1rnlx&e=1&subfolder_nav_tracking=1&dl=0>`_.

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
========

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

      *Analysis/Running options*
      --cutoff [number]     Minimum number of uncorrected UMIs required to process a barcode
                            DEFAULT: 1000
      --threads [number]    Number of processors/threads allocated.
                            DEFAULT:10

      *Other options*
      -h | --help           Show help message
      --version             Show program's version and exit


Example Command Lines
---------------------

::

    TEsingle  --threads 10 --stranded forward -b RNAseq.bam --GTF refseq_genes.gtf --TE rmsk_TE.gtf --project sample_test

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
--outSAMheaderHD @HD VN:1.4
--outSAMstrandField intronMotif
--outSAMtype BAM SortedByCoordinate
--sjdbScore 1
--soloType CB_samTagOut
--soloCBmatchWLtype 1MM
--soloCBwhitelist (10x whitelist)
--soloCellFilter Empty_Drops_CR 6000 0.99 10 45000 90000 500 0.01 20000 0.01 10000
--soloFeatures GeneFull


Copying & distribution
======================

TEsingle is part of `TEToolkit suite <https://www.mghlab.org/software>`_.

It is distributed under the BSD 3-clause license per ASAP Open Access (OA) policy, which facilitates the rapid and free exchange of scientific ideas and ensures that ASAP-funded research fund can be leveraged for future discoveries.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

A copy of BSD 3-clause licence is included along with the software, and can be accessed `here <https://github.com/mhammell-laboratory/TEsingle/blob/main/LICENSE>`_.


Acknowledgment
==============

- Contributors: Talitha Forcier, Oliver Tam, Cole Wunderlich & Molly Gale Hammell

This research was funded in part by Aligning Science Across Parkinson's (ASAP-000520) through the Michael J. Fox Foundation for Parkinson's Research (MJFF). Funding was also provided in part by the Chan-Zuckerberg Initiative (CZI) Neurodegeneration Challenge Network.
