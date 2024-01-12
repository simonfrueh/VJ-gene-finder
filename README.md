# VJ-gene-finder
## Extraction of putative chicken TCR VJ gene segments from a reference chromosome

VJ-gene-finder extracts functional V and J genes from the chicken genome. The algorithm is similar to the method used by Oliveri _et al._ [^1], but was modified to enable identification of chicken TCR V genes that are encoded by a single exon (together with the leader sequences) [^2] and an additional function was added to search and extract J gene candidates. The search parameters are based on characteristic biological patterns that define immunoglobulin V and J genes in many species and the nomenclature follows conventions of the international ImMunoGeneTics information system (IMGT)  [^3] 

The features used include:
1) Conserved amino acid residues at specific positions:
    - for V segments according to IMGT nomenclature: "1st-CYS", "CONSERVED-TRP" and “YYC/YFC/YLC/YHC/YIC/TFC” motif that includes the "2nd-CYS"
    - for J segments: “FG” motif
2) Conserved nucleic acid motifs in genes and at specific positions:
    - for V genes with a single-exon leader sequence: _ATG_ start codon
    - for V genes with a spliced leader sequence: splice acceptor sequence _AG_
    - for J segments: _TTYGGNNNNGG_ and _TNNBNRT_ and splice donor sequence _GTRDGD_
3) Conserved recombination signal sequences (RSS):
    - for both V and J segments: begin with _CAC_ nucleic acid motif

These motif requirements are combined with length constraints and the requirement for an open reading frame with or without splicing. 
Additionally, VJ-gene-finder tentatively assigns candidate V genes to chicken V gene families based on amino acid motifs near the 5’end: 

| TRV family name  | Amino acid motif |
| ---------------- | ---------------- |
| TRAV1            | QVQQ             |
| TRAV2            | VSQQ             |
| TRAV3            | LQYP             |
| TRBV1            | LQQT             |
| TRBV2            | EINQ             |
| TRBV3            | ITQW             |
| TRGV1            | QVLLQQ           |
| TRGV2            | PIQS             |
| TRGV3            | QAVPMQ or QAAPVQ |
| TRGV4            | LWQSP            |
| TRDV1            | ETSGGGV          |
| TRDV2            | LEASGGG          |
| TRDV3            | VEFGGDV          |
| TRDV4            | RIVEAG           |
| TRDV5            | EIHAKKSA         |
| TRDVH1           | QIEMVTT          |

Note: These motifs are based on chicken TRV families. If desired, the algorithm can be modified by altering the TRV | amino acid assignment in the v.py module. 

## Requirements
VJ-gene-finder requires that Python is installed and running.

## Installation
Download the latest release and unzip.

## Usage

```
python3 PATH_TO_VJ_GENE_FINDER_FOLDER/VJ_gene_finder.py PATH_TO_CHROMOSOME_FILE.fasta -o output
```

| Option                    | Argument                           | 
| ------------------------- | ---------------------------------- |
| filename (REQUIRED)       | provide input file in fasta format |
| -v / --version (optional) | use specific version               |
| -o / --output (optional)  | provide output directory           |

## Results

The VJ-gene-finder output generates three files. The files are named using the GenBank accession number that is part of the fasta header in files downloaded from the NCBI Genome assembly page.   

1) Accession_number.fasta
    - This file contains nucleic acid sequences of putative V and J genes with preliminary assignment to chicken TRV gene families and a unique identifier for each hit based on the location on the chromosome. If a TRV gene family could not be assigned, the hits are named either TRV with a unique ID, or TRV-SEL (single exon leader) with a unique ID when the V gene is part of an open reading frame with an _ATG_ start codon (within the length constraints). Additional information is provided in the fasta header, including the accession number, putative functionality (F = functional), and the start and end position on the forward or reverse strand (RC = reverse complement) of the original DNA sequence. 
2) Accession_number_aa.fasta
    - This file contains amino acid sequences of V and J gene hits with additional information in the fast header, as described above. 
3) Accession_number_RSS.fasta
    - This file contains recombination signal sequences (RSS) corresponding to V and J gene hits with additional information in the fasta header, as described above. 


**Results should be evaluated manually** because the search parameters are not specific to VJ genes only. 

1) Unspecific hits outside of the TCR loci can be removed.
    - The location of the TCR loci on the chromosome can be identified by alignment of the hit list, for example with [Clustal Omega] (https://www.ebi.ac.uk/Tools/msa/clustalo/).
    - Using a alignment viewer, such as [Jalview] (https://www.jalview.org/), the cluster of highly similar TCR sequences (named by gene family) on the forward and reverse strand can be distinguished from unspecific hits outside of the TCR loci.
2) The start and end position of V genes need to be evaluated.
    - Start position: If the V gene is part of an open reading frame with an _ATG_ start codon close to the putative V gene, VJ-gene-finder includes the whole open reading frame in the sequence. This behavior is desired for V genes with a single exon leader peptide, but provides false start positions for V genes with a spliced leader sequence. In this case the _AG_ splice site needs to be identified manually.  
    - End position: if there is a _CACAC_ motif, VJ-gene finder picks the first _CAC_ to define the end position, however, in our experience, the conserved RSS sequence usually begins at the second _CAC_.

[^1]: Olivieri, D., Faro, J., von Haeften, B., Sánchez-Espinel, C. and Gambón-Deza, F. (2013) ‘An automated algorithm for extracting functional immunologic V-genes from genomes in jawed vertebrates’, Immunogenetics, 65(9), pp. 691–702. Available at: https://doi.org/10.1007/s00251-013-0715-8.
[^2]: Göbel, T.W., Chen, C.L., Lahti, J., Kubota, T., Kuo, C.L., Aebersold, R., Hood, L. and Cooper, M.D. (1994) ‘Identification of T-cell receptor alpha-chain genes in the chicken’, Proceedings of the National Academy of Sciences of the United States of America, 91(3), pp. 1094–1098. Available at: https://doi.org/10.1073/pnas.91.3.1094.
[^3]: Lefranc, M.-P., Pommié, C., Ruiz, M., Giudicelli, V., Foulquier, E., Truong, L., Thouvenin-Contet, V. and Lefranc, G. (2003) ‘IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains’, Developmental & Comparative Immunology, 27(1), pp. 55–77. Available at: https://doi.org/10.1016/S0145-305X(02)00039-3.
