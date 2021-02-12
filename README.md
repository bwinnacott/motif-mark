# motif-mark

## Overview
This tool provides a visualization of motif sequences on gene sequences containing 
one exon and flanking introns. The current build can handle only gene sequences that 
have the exon denoted by capital nucleotides and the introns by lowercase nucleotides. 
All IUPAC degenerate base symbols (https://en.wikipedia.org/wiki/Nucleic_acid_notation) 
are considered in the motif sequences when locating matches in the gene sequence. 
Additionally, this tool can handle RNA motif sequences mapping to DNA gene sequences 
and vice versa. 

## Requirements
- Python 3.6+
- Pycairo 1.20.0

## Input File Requirements
- Sequence file: this program requires the input sequence file to be of type FASTA, 
with 1 exon (capitalized nucleotides) and flanking introns (lowercase)
- Motif file: a text file containing up to 5 motif sequences separated by line

## Arguments
The following arguments are required for running the program:
- ```-f```, ```--fa_file```: specifies absolute path to input FASTA file containing gene 
sequences
- ```-m```, ```--motif_file```: specifies absolute path to input motif sequence file
- ```-o```, ```--output_type```: type of image output desired; choices are 'png', 'svg', 
and 'pdf'

The following argument is optional:
- ```-c```, ```--colors```: specifies the color palette to use for visualizing different 
motif sequences on each gene; options are 'basic' (default; used when this argument is not 
called), 'earth', and 'colorblind_friendly'

## Example Output
![basic](https://github.com/bwinnacott/motif-mark/blob/main/example/Figure_1_earth.png)