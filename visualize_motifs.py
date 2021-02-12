#!/usr/bin/env python
import argparse
import cairo

def get_motifs(motif_file):
    '''
    Function to load the file containing the list of motif sequences to be searched in the 
    fasta file.
    motif_file: file containing list of motifs

    return: list of motifs
    '''
    # open connection to file
    motif = open(motif_file, 'r')
    # list for storing motif sequences
    motifs = []
    # iterate over each motif and add to list
    for line in motif:
        line = line.strip()
        motifs.append(line)
    # close the file
    motif.close()

    return motifs

def check_degen_symb(symbol):
    '''
    Function that takes in an IUPAC base symbol and returns the base(s) which it represents.
    symbol: IUPAC base symbol

    return: list of matching nts for input symbol
    '''
    # set the symbol mapping
    bases = {
        'A':['A','a'],'a':['A','a'],
        'C':['C','c'],'c':['C','c'],
        'G':['G','g'],'g':['G','g'],
        'T':['T','U','t','u'],'t':['T','U','t','u'],
        'U':['U','T','u','t'],'u':['U','T','u','t'],
        'W':['A','T','a','t'],'w':['A','T','a','t'],
        'S':['C','G','c','g'],'s':['C','G','c','g'],
        'M':['A','C','a','c'],'m':['A','C','a','c'],
        'K':['G','T','g','t'],'k':['G','T','g','t'],
        'R':['A','G','a','g'],'r':['A','G','a','g'],
        'Y':['C','T','c','t'],'y':['C','T','c','t'],
        'B':['C','G','T','c','g','t'],'b':['C','G','T','c','g','t'],
        'D':['A','G','T','a','g','t'],'d':['A','G','T','a','g','t'],
        'H':['A','C','T','a','c','t'],'h':['A','C','T','a','c','t'],
        'V':['A','C','G','a','c','g'],'v':['A','C','G','a','c','g'],
        'N':['A','C','G','T','a','c','g','t'],'n':['A','C','G','T','a','c','g','t'],
        'Z':[''],'z':['']
        }

    return bases[symbol]

def chars_to_skip(motif):
    '''
    This function takes in a motif pattern and generates an array containing the number of 
    characters to skip for each index position in the pattern, to optimize matching. The array 
    is returned. This is the preprocessing step for the KMP pattern matching algorithm.
    motif: sequence motif

    return: list of characters to skip for each letter in motif (when matching new window)
    '''
    # get length of motif
    l = len(motif)
    # initiate list of 0's of size l
    skip_chars = [0] * l
    # set counter to track index of current character in 'motif'; set another to track length 
    # of previous longest prefix that is also a suffix in 'motif' (i.e., matching sub-patterns)
    # Note: we start the first counter at 1, due to there not being a proper prefix for a single character
    i = 1
    prev_len = 0
    # loop until all characters are accounted for in 'motif'
    while i < l:
        # if characters match, add 1 to length counter and assign the current index in 'skip_chars' to 
        # the new length; increment the current index
        if motif[i] == motif[prev_len]:
            prev_len += 1
            skip_chars[i] = prev_len
            i += 1
        else:
            # if the longest prefix length is currently 0, assign current index in 'skip_chars' to 0, 
            # increment index counter
            if prev_len == 0:
                skip_chars[i] = 0
                i += 1
            # if previous longest length of prefix is not 0, assign new length (i.e., new pattern index for 
            # which to start comparing new window)
            else:
                prev_len = skip_chars[prev_len-1]

    return skip_chars

def search_sequence(motif,sequence):
    '''
    Function using KMP pattern matching algorithm to search a sequence for all instances of 
    current motif. All positions for which the motif is found are returned.
    motif: sequence motif
    sequence: sequence in fasta file

    return: list of indexes where a given motif matches the in the gene sequence
    '''
    # initiate list to store motif mapping positions
    motif_pos = []
    # get lengths of motif and sequence
    len_pat = len(motif)
    len_seq = len(sequence)
    # call 'chars_to_skip' function to get count of characters to be skipped
    skip_chars = chars_to_skip(motif)
    # initiate counters for sequence and motif index
    i = 0
    j = 0
    # loop until all characters are evaluated in sequence
    while i < len_seq:
        # if characters match, increment counters
        if sequence[i] in check_degen_symb(motif[j]):
            i += 1
            j += 1
        # if motif matches current window in sequence, get matching starting position; get 
        # index of next character in motif to be matched
        if j == len_pat:
            motif_pos.append(i-j)
            j = skip_chars[j-1]
        # if current characters don't match
        elif i < len_seq and sequence[i] not in check_degen_symb(motif[j]):
            # if first character of motif doesn't match window, shift window by 1
            if j == 0:
                i += 1
            # if current character (other than first) doesn't match, get index for next 
            # character in motif to be matched
            else:
                j = skip_chars[j-1]

    return motif_pos

def parse_fasta(fasta_file):
    '''
    Function that goes through fasta file and pulls out all relevant information for viewing 
    motifs on sequences.
    fasta_file: input fasta file containing gene sequences

    return: 
    num_seqs - number of sequences in fasta file
    longest_seq - longest sequence in fasta file
    sequences - dictionary (key: header, value: sequence)
    '''
    # initiate variables for figure dimension purposes
    num_seqs = 0
    longest_seq = 0
    # open connection to file
    with open(fasta_file,'r') as f:
        # initiate variables to store information
        header = None
        sequences = {}
        # iterate over each line
        for line in f:
            line = line.strip()
            # if header line
            if line.startswith('>'):
                num_seqs += 1
                # to account for None initially assigned to 'header'
                if header:
                    # get current longest sequence
                    longest_seq = max(longest_seq,len(sequences[header]))
                # store header
                sequences[line[1:]] = ''
                header = line[1:]
            # add the sequence line
            else:
                sequences[header] += line
        # account for the last sequence
        longest_seq = max(longest_seq,len(sequences[header]))
        
        return num_seqs,longest_seq,sequences

def get_color_palette(palette):
    '''
    Function used to get the RGB values for colors in the chosen palette.
    palette: string specified at command line for chosen color palette

    return: mapping of motif index to RGB color values
    '''
    # colorblind friendly palette
    if palette == 'colorblind_friendly':
        color_map = {0:(0.9,0.6,0),1:(0.35,0.7,0.9),2:(0,0.6,0.5),3:(0.95,0.9,0.25),4:(0,0.45,0.7)}
    # earth colors palette
    elif palette == 'earth':
        color_map = {0:(0.5,0.54,0.44),1:(0.81,0.7,0.53),2:(0.58,0.32,0.21),3:(0.16,0.27,0.25),4:(0.72,0.58,0.32)}
    # basic color palette (default)
    else:
        color_map = {0:(1,0,0),1:(0,1,0),2:(0,0,1),3:(0,1,1),4:(1,0,1)}

    return color_map

def draw_key(motifs,surface,color_map):
    '''
    Function to draw the figure legend at the top of the output image.
    motifs: list of motifs from input motif file
    surface: surface for which image is drawn
    color_map: dictionary containing rgb tuples specifying motif colors

    return: no return (modifies surface)
    '''
    # draw the introns and exon
    context = cairo.Context(surface)
    context.rectangle(15,25,30,30)
    context.fill()
    context.move_to(65,45)
    context.show_text('Exon')
    context.set_line_width(3)
    context.move_to(105,42)
    context.line_to(135,42)
    context.stroke()
    context.move_to(155,45)
    context.show_text('Intron')
    # variables to keep track of where motif labels should start in key 
    # and color assignment
    current_key_length = 200
    col_ind = 0
    # loop over motifs from longest to smallest
    for x in sorted(motifs,key=lambda x:len(x),reverse=True):
        # draw a colored rectangle for each
        context.rectangle(current_key_length,25,15,30)
        color = color_map[col_ind]
        (r,g,b) = color
        context.set_source_rgb(r,g,b)
        context.fill()
        # adjust position of text
        current_key_length += 25
        context.move_to(current_key_length,45)
        context.set_source_rgb(0,0,0)
        # add motif sequence
        context.show_text(x)
        # depending on length of motif, adjust position of next motif label
        if len(x) <= 2:
            current_key_length += 25
        elif 2<len(x)<=4:
            current_key_length += 40
        elif 4<len(x)<=6:
            current_key_length += 55
        elif 6<len(x)<=8:
            current_key_length += 70
        else:
            current_key_length += 85
        col_ind += 1
    # add a box around the legend
    context.rectangle(5,10,current_key_length-5,60)
    context.stroke()
    
def add_sequence_header(header,gene_ind,surface):
    '''
    Function to add modified header to top of each drawn sequence in fasta file.
    header: sequence header in fasta file
    gene_ind: index of gene in fasta file (1st sequence is 0, 2nd is 1, etc.)
    surface: surface for which image is drawn

    return: no return (modifies surface)
    '''
    # set position on surface for header
    context = cairo.Context(surface)
    context.set_source_rgb(0,0,0)
    context.move_to(5,gene_ind*100+10)
    # modify the header to make it more readable
    header = header.split(' ')
    context.show_text('Gene: ' + str(header[0]) + 2*' ' + '-->' + 2*' ' + 
                      'Chromosome: ' + header[1].split(':')[0] + 2*' ' + 
                      '-->' + 2*' ' + 'Coordinates: ' + str(header[1].split(':')[1]))
    
def add_introns_exon(sequence,gene_ind,surface):
    '''
    Function to add the introns and exon for a given sequence.
    sequence: gene sequence
    gene_ind: index of current gene sequence in fasta file
    surface: surface for which image is drawn

    return: no return (modifies surface)
    '''
    context = cairo.Context(surface)
    # get the exon position
    exon_pos = [i for i,char in enumerate(sequence) if char.isupper()]
    context.set_line_width(3)
    # draw the intron on the surface
    context.move_to(5,gene_ind*100+50)
    context.line_to(len(sequence)+5,gene_ind*100+50)
    context.set_source_rgb(0,0,0)
    context.stroke()
    # draw the exon on the surface
    context.rectangle(min(exon_pos)+5,25+100*gene_ind,max(exon_pos)-min(exon_pos)+1,50)
    context.set_source_rgb(0,0,0)
    context.fill()
    
def get_motif_pos(motifs,seq):
    '''
    Function to obtain the positions on the sequence where motifs match.
    motifs: list of motifs
    seq: sequence from fasta file

    return: dictionary of motif matching indexes
    '''
    # initiate dictionary
    motif_mapping = {}
    # loop over each motif and add motif as key and list of indexes in sequence where 
    # it matches as value
    for motif in motifs:
        motif_mapping[motif] = search_sequence(motif,seq)
        
    return motif_mapping

def draw_motif(motif,motif_map,color_map,gene_ind,motif_ind,num_motifs,surface):
    '''
    Function used to draw the motif sequences on the genes. Motifs are offset to allow 
    for better visualization of overlap.
    motif: motif sequence
    motif_map: dictionary mapping motif sequences to index positions on gene sequence
    color_map: specifies color mapping for motifs
    gene_ind: index for gene in fasta file
    motif_ind: index for motif
    num_motifs: number of motifs in motif input file
    surface: surface for drawing image

    return: no return (modifies object)
    '''
    # get the color for the current motif
    context = cairo.Context(surface)
    color = color_map[motif_ind]
    (r,g,b) = color
    # get factors for adjusting the position of the motif on the image
    pos_factor = 0 if num_motifs == 1 else 25/(num_motifs-1)
    len_factor = 50 if num_motifs == 1 else 25
    # draw the motif for each location it matched to the gene sequence
    for pos in motif_map[motif]:
        context.rectangle(pos+5,25+(pos_factor*motif_ind)+100*gene_ind,len(motif),len_factor)
        context.set_source_rgba(r,g,b)
        context.fill()
        
def output_image(input_filename,width,height,out_filetype):
    '''
    Function to write out the image to a specific file type.
    input_filename: input fasta file name
    width: width of output image
    height: height of output image
    out_filetype: file type specified at command line

    return: surface for drawing image
    '''
    # svg output
    if out_filetype == 'svg':
        return cairo.SVGSurface(str(input_filename.split('.')[0] + '.svg'),width,height)
    # png output
    elif out_filetype == 'png':
        return cairo.ImageSurface(cairo.FORMAT_ARGB32,width,height)
    # pdf output
    else:
        return cairo.PDFSurface(str(input_filename.split('.')[0] + '.pdf'),width,height)

def main():
    # set command line arguments
    parser = argparse.ArgumentParser(description='Tool used to mark motifs on gene sequences containing 1 exon and introns flanking either side. \
                                    Can handle a maximum of 5 motifs on each sequence.''')
    parser.add_argument('-f','--fa_file',action='store',required=True,type=str,help='Specifies the input FASTA file for which sequence motif marking is desired.')
    parser.add_argument('-m','--motif_file',action='store',required=True,type=str,help='Specifies the input .txt file containing the list of motifs.')
    parser.add_argument('-o','--output_type',action='store',required=True,choices=['svg','pdf','png'],type=str,help='Designates the output file type.')
    parser.add_argument('-c','--colors',action='store',choices=['earth','colorblind_friendly'],type=str,help='Specifies color palette used for differentiating motifs. \
                        If not provided at command line, the "basic" color palette is used as default.')
    # extract the argument information
    args = parser.parse_args()
    # get fasta info and motif list
    num_seq,longest_seq,sequences = parse_fasta(args.fa_file)
    motifs = get_motifs(args.motif_file)
    # get color palette and image for drawing
    color_map = get_color_palette(args.colors)
    width,height = longest_seq+10,100*(num_seq+1)
    surface = output_image(args.fa_file,width,height,args.output_type)
    # draw the key
    draw_key(motifs,surface,color_map)
    gene_ind = 1
    # draw the sequences one at a time
    for head,seq in sequences.items():
        add_sequence_header(head,gene_ind,surface)
        add_introns_exon(seq,gene_ind,surface)
        motif_mapping = get_motif_pos(motifs,seq)
        motif_ind = 0
        # get all motifs marked for current sequence
        for motif in sorted(motif_mapping,key=lambda x:len(x),reverse=True):
            draw_motif(motif,motif_mapping,color_map,gene_ind,motif_ind,len(motifs),surface)
            motif_ind += 1
        gene_ind += 1
    # finishing the surface object and writing out to file type specified
    if args.output_type == 'png':
        surface.write_to_png(str('Figure_1.fasta'.split('.')[0] + '.png'))
    else:
        surface.finish()

if __name__ == "__main__":
    main()