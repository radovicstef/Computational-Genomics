# This is a Pair-end sequencer simulator, project for Computer Genomics class


# Methods
import numpy.random
from random import randrange


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
    fastq_1 = open("genome_1.fastq", "w")
    fastq_2 = open("genome_2.fastq", "w")
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
            read_2 = complement_read(read_2)
            read_2 = reverse_read(read_2)
            read_2_quality = get_quality(avg_quality, SIGMA, read_length)
            read_2_quality = reverse_read(read_2_quality)

            fastq_1.write("{}\n{}\n+\n{}\n".format(read_1_id, read_1, read_1_quality))
            fastq_2.write("{}\n{}\n+\n{}\n".format(read_2_id, read_2, read_2_quality))


# Global variables
genome = {}  # dictionary - sequence_name:sequence from FASTA file
genome_length = 0
num_of_reads = 0  # number of reads = number of fragments
num_of_pair_end_reads = 0  # number of pair-end reads = 2 * number of reads
SIGMA = 3


# The main program
def simulate(genome_file, avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv=0, prob_ins=0, prob_del=0):
    global genome
    read_genome(genome_file)
    print("Genome size from FASTA file: {}".format(genome_length))


if __name__ == "__main__":
    simulate("genomeSample.fa", 100, 4, 7, 10)
