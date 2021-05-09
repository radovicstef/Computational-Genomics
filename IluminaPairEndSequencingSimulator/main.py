# This is a Pair-end sequencer simulator, project for Computer Genomics class


# Methods
import numpy as np
import numpy.random
from random import randrange
import re


def check_input_parameters(avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv, prob_ins, prob_del):
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
    global genome, genome_length
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
    genome_length = len(genome)


def calculate_number_of_reads(sequence_length, coverage, read_length):
    global num_of_reads, num_of_pair_end_reads
    # Number of reads = Number of fragments
    num_of_reads = int(round(coverage * sequence_length / (read_length * 2)))  # Lander-Waterman formula
    num_of_pair_end_reads = 2 * num_of_reads


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
    for sequence_name in genome:
        sequence_length = len(genome[sequence_name])
        snv_num = round(sequence_length * prob_snv)
        ins_num = round(sequence_length * prob_ins)
        del_num = round(sequence_length * prob_del)

        for i in range(snv_num):
            mutation_position = randrange(sequence_length-1)
            genome[sequence_name] = genome[sequence_name][0:mutation_position] + get_mutated_single_nucleotide(genome[sequence_name][mutation_position]) + genome[sequence_name][mutation_position+1:]

        for i in range(ins_num):
            mutation_position = randrange(sequence_length - 1)
            genome[sequence_name] = genome[sequence_name][0:mutation_position] + bases[randrange(4)] + genome[sequence_name][mutation_position:]

        for i in range(del_num):
            mutation_position = randrange(sequence_length - 1)
            genome[sequence_name] = genome[sequence_name][0:mutation_position] + genome[sequence_name][mutation_position+1:]


def get_reads(avg_quality, coverage, read_length, insert_size):
    global SIGMA
    global fragment_positions
    fastq_1 = open("genome_1.fastq", "w")
    fastq_2 = open("genome_2.fastq", "w")
    sam_file = open("final_sam_file.sam", "w")
    for sequence_name, sequence in genome.items():
        calculate_number_of_reads(len(sequence), coverage, read_length)
        # Fragment position should be uniformly distributed to suit specified coverage
        fragment_positions = numpy.random.uniform(0, len(sequence) - insert_size, num_of_reads)
        for i in range(num_of_reads):
            fragment_position = round(fragment_positions[i])
            if fragment_position == (len(sequence) - insert_size):
                fragment_position -= 1

            read_1_id = "@" + sequence_name + str(i) + "/1"
            read_1 = sequence[fragment_position:(fragment_position + read_length)]
            read_1_quality = get_quality(avg_quality, SIGMA, read_length)

            read_2_id = "@" + sequence_name + str(i) + "/2"
            read_2 = sequence[(fragment_position + insert_size - read_length):(fragment_position + insert_size)]
            read_2_sam = read_2
            read_2 = complement_read(read_2)
            read_2 = reverse_read(read_2)
            read_2_quality = get_quality(avg_quality, SIGMA, read_length)
            read_2_quality_sam = read_2_quality
            read_2_quality = reverse_read(read_2_quality)

            fastq_1.write("{}\n{}\n+\n{}\n".format(read_1_id, read_1, read_1_quality))
            fastq_2.write("{}\n{}\n+\n{}\n".format(read_2_id, read_2, read_2_quality))
            sam_file.write("{}\t{}\t{}\t{}\n".format(read_1_id, fragment_position + 1, read_1, read_1_quality)) #+1 because it's 1-based
            sam_file.write("{}\t{}\t{}\t{}\n".format(read_2_id, fragment_position + insert_size - read_length + 1, read_2_sam, read_2_quality_sam))


# Compare BWA-MEM with generated SAM file
def compare_sam_bwa_mem(bwa_mem_sam_path, generated_sam_path):
    all_reads = 0
    matched = 0
    skip_lines = 0
    with open(bwa_mem_sam_path, "r") as bwa_mem_sam:
        for line in bwa_mem_sam:
            if line[0] == '@':
                skip_lines += 1
                continue
            else:
                bwa_mem_sam.close()
                break
    with open(bwa_mem_sam_path, "r") as bwa_mem_sam, open(generated_sam_path, "r") as generated_sam:
        for i in range(skip_lines):
            line = bwa_mem_sam.readline()
        for line1, line2 in zip(bwa_mem_sam, generated_sam):
            line1_split = re.split(r'\t+', line1)
            line2_split = re.split(r'\t+', line2)
            all_reads += 1
            if line1_split[3] == line2_split[1]:
                matched = matched + 1
    print("BWA-MEM efficiency: " + str(matched/all_reads))
    return matched/all_reads


def compare_sam_bowtie(bowtie_sam_path, generated_sam_path):
    all_reads = 0
    matched = 0
    bowtie_reads = {}
    with open(bowtie_sam_path, "r") as bowtie_sam:
        for line in bowtie_sam:
            if line[0] == '@':
                continue
            line_split = re.split(r'\t+', line)
            bowtie_reads['@' + line_split[0]] = line_split[3]
    with open(generated_sam_path, "r") as generated_sam:
        for line in generated_sam:
            line_split = re.split(r'\t+', line)
            all_reads += 1
            if bowtie_reads[line_split[0]] == line_split[1]:
                matched += 1
    print("Bowtie efficiency: " + str(matched/all_reads))
    return matched/all_reads


# Global variables
genome = {}  # dictionary - sequence_name:sequence from FASTA file
genome_length = 0
num_of_reads = 0  # number of reads = number of fragments
num_of_pair_end_reads = 0  # number of pair-end reads = 2 * number of reads
SIGMA = 3
fragment_positions = []


# The main program
def simulate(genome_file, avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv=0, prob_ins=0, prob_del=0):
    global genome
    valid = check_input_parameters(avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv, prob_ins, prob_del)
    if valid == 0:
        print("The input parameters are not valid!")
        return
    read_genome(genome_file)
    #get_reads(30, 10, 10, 25)
    #compare_sam_bwa_mem("genome.sam", "final_sam_file.sam")
    compare_sam_bowtie("genome_bowtie.sam", "final_sam_file.sam")
    #nesto


if __name__ == "__main__":
    simulate("genomeSample.fa", 30, 10, 10, 25)
