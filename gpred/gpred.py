"""Gene prediction from a bacterial genome based on Shine-Dalgarno motifs."""

import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__,
                                 usage=f"{sys.argv[0]} -h")
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path,
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome.
    """
    sequence = ""
    with open(fasta_file, "r", encoding="utf-8") as fich:
        for line in fich:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None.
    """
    match = start_regex.search(sequence, start, stop)
    if match:
        return match.start()
    return None


def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None.
    """
    for i in range(start, len(sequence) - 2, 3):
        if stop_regex.match(sequence, i):
            return i
    return None


def has_shine_dalgarno(shine_regex: Pattern, sequence: str, start: int,
                        max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of
    the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    debut = start - max_shine_dalgarno_distance
    if debut < 0:
        return False
    return shine_regex.search(sequence, debut, start - 6) is not None


def predict_genes(sequence: str, start_regex: Pattern, stop_regex: Pattern, shine_regex: Pattern,
                  min_gene_len: int, max_shine_dalgarno_distance: int, min_gap: int) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of
    the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    probable_genes = []
    current_pos = 0
    while len(sequence) - current_pos >= min_gap:
        start = find_start(start_regex, sequence, current_pos, len(sequence))
        if start is None:
            break
        stop = find_stop(stop_regex, sequence, start)
        if stop is not None:
            gene_len = stop + 3 - start
            if gene_len >= min_gene_len:
                if has_shine_dalgarno(shine_regex, sequence, start,
                                      max_shine_dalgarno_distance):
                    probable_genes.append([start + 1, stop + 3])
                    current_pos = stop + 3 + min_gap
                    continue
        current_pos = start + 1
    return probable_genes


def write_genes_pos(predicted_genes_file: Path, probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as csv_file:
            writer = csv.writer(csv_file, delimiter=",")
            writer.writerow(["Start", "Stop"])
            writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def write_genes(fasta_file: Path, sequence: str, probable_genes: List[List[int]], sequence_rc: str,
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of
    each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt", encoding="utf-8") as fasta:
            nb_genes = len(probable_genes)
            for i, gene_pos in enumerate(probable_genes):
                fasta.write(f">gene_{i+1}{os.linesep}"
                            f"{textwrap.fill(sequence[gene_pos[0]-1:gene_pos[1]])}{os.linesep}")
            for j, gene_pos in enumerate(probable_genes_comp):
                fasta.write(f">gene_{nb_genes+1+j}{os.linesep}"
                            f"{textwrap.fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])}{os.linesep}")
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


#==============================================================
# Main program
#==============================================================
def main() -> None: # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()
    # Let us do magic in 5' to 3'
    sequence = read_fasta(args.genome_file)
    probable_genes = predict_genes(sequence, start_regex, stop_regex, shine_regex,
                                   args.min_gene_len,
                                   args.max_shine_dalgarno_distance,
                                   args.min_gap)

    # Don't forget to uncomment !!!
    # Call these function in the order that you want
    # We reverse and complement
    sequence_rc = reverse_complement(sequence)
    probable_genes_comp = predict_genes(sequence_rc, start_regex, stop_regex,
                                        shine_regex, args.min_gene_len,
                                        args.max_shine_dalgarno_distance,
                                        args.min_gap)
    # Call to output functions
    write_genes_pos(args.predicted_genes_file, probable_genes)
    write_genes(args.fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp)



if __name__ == '__main__':
    main()
