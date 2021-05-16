# This is a Pair-end sequencer simulator, project for Computer Genomics class


# Methods
import numpy as np
import numpy.random
from random import randrange
import re
import argparse


def check_input_parameters(avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv=0, prob_ins=0, prob_del=0):
    if avg_nucleotide_quality < 33 or avg_nucleotide_quality > 73:
        return 0
    if read_length < 0 or insert_size < 0:
        return 0
    if prob_snv < 0 or prob_snv > 1:
        return 0
    if prob_ins < 0 or prob_ins > 1:
        return 0
    if prob_del < 0 or prob_del > 1:
        return 0
    if read_length > insert_size:
        return 0
    return 1


def read_genome(file_name): 
    global genome
    print("Reading reference genome...")
    genome_file = open(file_name, "r")
    genome_name = ""
    lines = genome_file.read().splitlines()  # returns the list of lines without \n
    for line in lines:
        # get the name of the sequence
        if line.find('>') != -1:
            genome_name = line.split()[0][1:]
            genome[genome_name] = ""
        # get the sequence
        else:
            genome[genome_name] = genome[genome_name] + line

    genome_file.close()


def calculate_number_of_reads(sequence_length, coverage, read_length):
    # Number of reads = Number of fragments
    num_of_reads = int(round(coverage * sequence_length / (read_length * 2)))  # Lander-Waterman formula
    num_of_pair_end_reads = 2 * num_of_reads
    return num_of_reads, num_of_pair_end_reads


# Second read should be reversed, as well as the quality string
def reverse_read(read):
    return read[::-1]


# Second read should be complemented
def complement_read(read):
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N'
    }
    return ''.join([complement[base] for base in read])


# Quality range = 0 - 40, ASCII characters range = 33 - 73
# Illumina reference: https://bit.ly/3eb8xw0
def get_quality(avg_quality, sigma, read_length):
    qualities = [round(quality) for quality in numpy.random.normal(avg_quality, sigma, read_length)]
    for quality in qualities:
        if quality < 0:
            quality = 0
        if quality > 40:
            quality = 40
    return "".join([chr(quality+33) for quality in qualities])


def get_mutated_single_nucleotide(nucleotide):
    bases = ["A", "C", "G", "T", "N"]
    bases.remove(nucleotide)
    return bases[randrange(4)]


def add_mutations(prob_snv=0, prob_ins=0, prob_del=0):
    bases = ['A', 'C', 'G', 'T']
    print("Adding mutations...")
    for sequence_name in genome:
        sequence_length = len(genome[sequence_name])
        snv_num = round(sequence_length * prob_snv)
        ins_num = round(sequence_length * prob_ins)
        del_num = round(sequence_length * prob_del)

        for i in range(snv_num):
            mutation_position = randrange(sequence_length-1)
            genome[sequence_name] = genome[sequence_name][0:mutation_position] + get_mutated_single_nucleotide(genome[sequence_name][mutation_position]) + genome[sequence_name][mutation_position+1:]

        for i in range(ins_num):
            mutation_position = randrange(len(genome[sequence_name]) - 1)
            genome[sequence_name] = genome[sequence_name][0:mutation_position] + bases[randrange(4)] + genome[sequence_name][mutation_position:]

        for i in range(del_num):
            mutation_position = randrange(len(genome[sequence_name]) - 1)
            genome[sequence_name] = genome[sequence_name][0:mutation_position] + genome[sequence_name][mutation_position+1:]


# Generate reads, output FASTQ files, SAM file
def get_reads_and_generate_files(avg_quality, coverage, read_length, insert_size):
    global SIGMA
    fastq_1 = open("genome_1.fastq", "w")
    fastq_2 = open("genome_2.fastq", "w")
    sam_file = open("final_sam_file.sam", "w")
    print("Sequencing and generating FASTQ, SAM files...")
    for sequence_name, sequence in genome.items():
        num_of_reads, num_of_pair_end_reads = calculate_number_of_reads(len(sequence), coverage, read_length)
        # Fragment position should be uniformly distributed to suit specified coverage
        fragment_positions = numpy.random.uniform(0, len(sequence) - insert_size, num_of_reads)
        for i in range(num_of_reads):
            fragment_position = round(fragment_positions[i])
            if fragment_position == (len(sequence) - insert_size):  # Edge case
                fragment_position -= 1

            read_1_id = "@" + sequence_name + str(i) + "/1"
            read_1 = sequence[fragment_position:(fragment_position + read_length)]
            read_1_quality = get_quality(avg_quality, SIGMA, read_length)

            read_2_id = "@" + sequence_name + str(i) + "/2"
            read_2 = sequence[(fragment_position + insert_size - read_length):(fragment_position + insert_size)]
            read_2_sam = read_2  # Read2 shouldn't be reversed complemented in SAM file
            read_2 = complement_read(read_2)
            read_2 = reverse_read(read_2)
            read_2_quality = get_quality(avg_quality, SIGMA, read_length)
            read_2_quality_sam = read_2_quality  # Read2 quality shouldn't be reversed in SAM file
            read_2_quality = reverse_read(read_2_quality)

            fastq_1.write("{}\n{}\n+\n{}\n".format(read_1_id, read_1, read_1_quality))
            fastq_2.write("{}\n{}\n+\n{}\n".format(read_2_id, read_2, read_2_quality))
            sam_file.write("{}\t{}\t{}\t{}\n".format(read_1_id, fragment_position + 1, read_1, read_1_quality))  # +1 because SAM file is 1-based
            sam_file.write("{}\t{}\t{}\t{}\n".format(read_2_id, fragment_position + insert_size - read_length + 1, read_2_sam, read_2_quality_sam))


# Compare BWA-MEM with generated SAM file
def compare_sam(tool_sam_path, generated_sam_path):
    all_reads = 0
    matched = 0
    skip_lines = 0
    print("Comparing SAM files...")
    with open(tool_sam_path, "r") as tool_sam:
        for line in tool_sam:
            if line[0] == '@':  # Skip header
                skip_lines += 1
                continue
            else:
                break
    with open(tool_sam_path, "r") as tool_sam, open(generated_sam_path, "r") as generated_sam:
        for i in range(skip_lines):
            line = tool_sam.readline()
        for line1, line2 in zip(tool_sam, generated_sam):
            line1_split = re.split(r'\t+', line1)
            line2_split = re.split(r'\t+', line2)
            all_reads += 1
            if line1_split[3] == line2_split[1]:
                matched = matched + 1
            else:
                print("Unmatched")
                print("-----")
                print(line1)
                print(line2)
                print("-----")
    print("BWA-MEM efficiency: " + str(matched/all_reads))
    print("Finished!")
    return matched/all_reads


# Global variables
genome = {}  # dictionary - sequence_name:sequence from FASTA file, used in test.py and therefore global
SIGMA = 3  # constant value, therefore global


# The main program
def simulate(genome_file, avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv=0, prob_ins=0, prob_del=0):
    global genome
    print("Starting simulation...")
    if None in (genome_file, avg_nucleotide_quality, coverage, read_length, insert_size):
        print("Please, make sure you specified all necessary parameters - refgenome, avgquality, coverage, readlength, insertsize")
        return
    valid = check_input_parameters(avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv, prob_ins, prob_del)
    if valid == 0:
        print("The input parameters are not valid!")
        return
    read_genome(genome_file)
    add_mutations(prob_snv, prob_ins, prob_del)
    get_reads_and_generate_files(avg_nucleotide_quality, coverage, read_length, insert_size)
    print("Finished!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Illumina Pair-End Sequencing Simulator"
    )  # Initializing parser

    # Add the parameters positional/optional
    parser.add_argument('func', help="Choose functionality: 1 - generate FASTQ files and SAM file, 2 - Compare BWA-MEM, 3 - Compare Bowtie", type=int)
    parser.add_argument('--refgenome', help="Path to reference genome", type=str)
    parser.add_argument('--avgquality', help="Average nucletiode quality", type=int)
    parser.add_argument('--coverage', help="Coverage", type=int)
    parser.add_argument('--readlength', help="Read length", type=int)
    parser.add_argument('--insertsize', help="Insert size", type=int)
    parser.add_argument('--probsnv', help="Probability snv mutation", type=float, default=0)
    parser.add_argument('--probins', help="Probability insert mutation", type=float, default=0)
    parser.add_argument('--probdel', help="Probability delete mutation", type=float, default=0)
    parser.add_argument('--bwapath', help="Path to BWA-MEM tool generated SAM file", type=str)
    parser.add_argument('--bowtiepath', help="Path to Bowtie tool generated SAM file", type=str)
    parser.add_argument('--sampath', help="Path to simulator generated SAM file", type=str)

    args = parser.parse_args()  # Parse the arguments

    if args.func == 1:  # Simulate Illumina pair-end sequencing
        simulate(args.refgenome, args.avgquality, args.coverage, args.readlength, args.insertsize, args.probsnv, args.probins, args.probdel)
    if args.func == 2:  # Compare generated SAM file with BWA-MEM tool generated SAM file
        if None in (args.bwapath, args.sampath):
            print("Please, make sure you specified all necessary parameters - bwapath, sampath")
        else:
            compare_sam(args.bwapath, args.sampath)
    if args.func == 3:  # Compare generated SAM file with Bowtie tool generated SAM file
        if None in (args.bowtiepath, args.sampath):
            print("Please, make sure you specified all necessary parameters - bowtiepath, sampath")
        else:
            compare_sam(args.bowtiepath, args.sampath)
