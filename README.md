# *VJ-gene-finder*

## Extraction of putative chicken TCR VJ gene segments from a reference chromosome

*VJ-gene-finder* extracts functional V and J genes from the chicken genome. The algorithm is similar to the method used by Oliveri *et al.* <sup>[\[1\]](#fn1)</sup>, but was modified to enable identification of chicken TCR V genes that are encoded by a single exon (together with the leader sequences) <sup>[\[2\]](#fn2)</sup> and an additional function was added to search and extract J gene candidates. The search parameters are based on characteristic biological patterns that define immunoglobulin V and J genes in many species and the nomenclature follows conventions of the international ImMunoGeneTics information system (IMGT) <sup>[\[3\]](#fn3)</sup>

The features used include:

1.  Conserved amino acid residues at specific positions:
    - for V segments according to IMGT nomenclature: "1st-CYS", "CONSERVED-TRP" and “YYC/YFC/YLC/YHC/YIC/TFC” motif that includes the "2nd-CYS"
    - for J segments: “FG” motif
2.  Conserved nucleic acid motifs in genes and at specific positions:
    - for V genes with a single exon leader sequence: *ATG* start codon
    - for V genes with a spliced leader sequence: splice acceptor sequence *AG*
    - for J segments: *TTYGGNNNNGG* and *TNNBNRT* and splice donor sequence *GTRDGD*
3.  Conserved recombination signal sequences (RSS):
    - for both V and J segments: begin with *CAC* nucleic acid motif

These motif requirements are combined with length constraints and the requirement for an open reading frame. For more information please refer to the *VJ-gene-finder* publication.

Additionally, *VJ-gene-finder* tentatively assigns candidate V genes to chicken V gene families based on amino acid motifs near the 5’end.  
**Note:** These motifs were derived from chicken TRV families. If desired, the algorithm can be modified by altering the TRV | amino acid assignment in the v.py module. Alternatively, the assignment can be turned off.

| TRV family name | Amino acid motif |
| --- | --- |
| TRAV1 | QVQQ |
| TRAV2 | VSQQ |
| TRAV3 | LQYP |
| TRBV1 | LQQT |
| TRBV2 | EINQ |
| TRBV3 | ITQW |
| TRGV1 | QVLLQQ |
| TRGV2 | PIQS |
| TRGV3 | QAVPMQ or QAAPVQ |
| TRGV4 | LWQSP |
| TRDV1 | ETSGGGV |
| TRDV2 | LEASGGG |
| TRDV3 | VEFGGDV |
| TRDV4 | RIVEAG |
| TRDV5 | EIHAKKSA |
| TRDVH1 | QIEMVTT |

## Requirements

*VJ-gene-finder* requires that Python and [Biopython](https://biopython.org/) are installed and running.

## Installation

Download the latest release and unzip.

## Usage

```
python3 VJ_gene_finder.py [-h] [-v] [-o OUTPUT] [-st] [-ss] filename
```

| Argument |     |
| --- | --- |
| filename | Input file in fasta format |

| Options |     |
| --- | --- |
| \-h, --help | Show help message and exit |
| \-v, --version | Show version number and exit |
| \-o OUTPUT, --output OUTPUT | Output directory |
| \-st, --skip_trgf | Skip assignment of TR family name (based on amino acid motif) for V genes  <br>and use default names instead  <br>(TRV SEL for V genes with single-exon leader peptide,  <br>TRV for V genes with two-exon leader peptide). |
| \-ss, --skip_selp | Skip search for V genes with single-exon leader peptide. |

## Results

The *VJ-gene-finder* output generates three files. The files are named using the GenBank accession number that is part of the fasta header in files downloaded from the NCBI Genome assembly page.

1.  Accession_number.fasta
    - This file contains nucleic acid sequences of putative V and J genes with preliminary assignment to chicken TRV gene families (optional) and a unique identifier for each hit based on the position in the query sequence. If a TRV gene family could not be assigned, the hits are named either TRV with a unique ID, or TRV-SEL (single exon leader) with a unique ID when the V gene is part of an open reading frame with an *ATG* start codon (within the length constraints). Additional information is provided in the fasta header, including the accession number, putative functionality (F = functional), and the start and end position on the forward or reverse strand (RC = reverse complement) of the query DNA sequence.
2.  Accession_number_aa.fasta
    - This file contains amino acid sequences of V gene hits with additional information in the fasta header, as described above.
3.  Accession_number_RSS.fasta
    - This file contains recombination signal sequences (RSS) corresponding to V and J gene hits with additional information in the fasta header, as described above.

**Results should be evaluated manually** because the search parameters are not specific to VJ genes only.

1.  Unspecific hits outside of the TCR loci can be removed.
    - The location of the TCR loci on the chromosome can be identified by alignment of the hit list, for example with [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/).
    - Using a alignment viewer, such as [Jalview](https://www.jalview.org/), the cluster of highly similar TCR sequences (named by gene family) on the forward and reverse strand can be distinguished from unspecific hits outside of the TCR loci.
2.  The start and end position of V genes need to be evaluated.
    - Start position: If the V gene is part of an open reading frame with an *ATG* start codon close to the putative V gene, VJ-gene-finder includes the whole open reading frame in the sequence. This behavior is desired for V genes with a single exon leader peptide, but provides false start positions for V genes with a spliced leader sequence. In this case the *AG* splice site needs to be identified manually.
    - End position: if there is a *CACAC* motif, VJ-gene finder picks the first *CAC* to define the end position, however, in our experience, the conserved RSS sequence usually begins at the second *CAC*.

&nbsp;  
For more information please refer to the *VJ-gene-finder* publication.

## Query sequences from other species

*VJ-gene-finder* relies on conserved biological features that are not chicken-specific per se, so it is likely that the algorithm would also identify TCR V genes from other species. Nonetheless, the search parameters may need to be adjusted for query sequences from non-chicken species. For additional details about the algorithm and features that may be modified, please refer to the *VJ-gene-finder* publication. For species without single exon leader peptides, this search feature can be turned off. In addition, the assignment to TRV families can either be modified, as described above, or disabled entirely for species where information about V families is not available.

* * *

1.  Olivieri, D., Faro, J., von Haeften, B., Sánchez-Espinel, C. and Gambón-Deza, F. (2013) ‘An automated algorithm for extracting functional immunologic V-genes from genomes in jawed vertebrates’, Immunogenetics, 65(9), pp. 691–702. Available at: https://doi.org/10.1007/s00251-013-0715-8. [↩︎](#fnref1)
    
2.  Göbel, T.W., Chen, C.L., Lahti, J., Kubota, T., Kuo, C.L., Aebersold, R., Hood, L. and Cooper, M.D. (1994) ‘Identification of T-cell receptor alpha-chain genes in the chicken’, Proceedings of the National Academy of Sciences of the United States of America, 91(3), pp. 1094–1098. Available at: https://doi.org/10.1073/pnas.91.3.1094. [↩︎](#fnref2)
    
3.  Lefranc, M.-P., Pommié, C., Ruiz, M., Giudicelli, V., Foulquier, E., Truong, L., Thouvenin-Contet, V. and Lefranc, G. (2003) ‘IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains’, Developmental & Comparative Immunology, 27(1), pp. 55–77. Available at: [https://doi.org/10.1016/S0145-305X(02)00039-3](https://doi.org/10.1016/S0145-305X%2802%2900039-3 "https://doi.org/10.1016/S0145-305X(02)00039-3"). [↩︎](#fnref3)
