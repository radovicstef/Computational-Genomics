import main

# Run main method
main.main()

# Test

# Test if sample genome is correctly extracted from FASTA file
valid_genome_file = open('genome.txt', 'r')
valid_genome = valid_genome_file.readline()
assert(valid_genome == main.genome)