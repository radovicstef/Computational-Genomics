# This is a Pair-end sequencer simulator, project for Computer Genomics class


# Methods
def read_genome(file_name):
    global genome
    genome_file = open(file_name, "r")
    lines = genome_file.read().splitlines() #returns the list of lines without \n
    for line in lines:
        # get only the nucleotides
        if line.find('>') == -1:
            genome = genome + line
    genome_file.close()


# Global variables
genome = ""


# The main program
def main():
    global genome
    read_genome("genomeSample.fa")
    print("Genome from FASTA file: {}".format(genome))

if __name__ == "__main__":
    main()
