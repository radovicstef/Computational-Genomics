# This is a Pair-end sequencer simulator, project for Computer Genomics class


# Methods
def read_genome(file_name):
    global genome, genome_length
    genome_file = open(file_name, "r")
    lines = genome_file.read().splitlines()  # returns the list of lines without \n
    for line in lines:
        # get only the nucleotides
        if line.find('>') == -1:
            genome = genome + line
    genome_file.close()
    genome_length = len(genome)


def calculate_number_of_reads(coverage, read_length):
    global genome_length, num_of_reads, num_of_pair_end_reads
    num_of_reads = int(round(coverage * genome_length / (read_length * 2)))  # Lander-Waterman formula
    num_of_pair_end_reads = 2 * num_of_reads


# Global variables
genome = ""
genome_length = 0
num_of_reads = 0  # number of reads = number of fragments
num_of_pair_end_reads = 0  # number of pair-end reads = 2 * number of reads


# The main program
def simulate(genome_file, avg_nucleotide_quality, coverage, read_length, insert_size, prob_snv=0, prob_ins_del=0):
    global genome
    read_genome(genome_file)
    calculate_number_of_reads(coverage, read_length)
    print("Number of reads: {}".format(num_of_reads))
    print("Genome size from FASTA file: {}".format(genome_length))


if __name__ == "__main__":
    simulate("genomeSample.fa", 100, 4, 7, 10)
