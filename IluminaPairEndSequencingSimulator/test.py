import main

# Run main method
main.simulate("genomeSample.fa", 100, 4, 7, 10)

# Tests

# Test if sample genome is correctly extracted from FASTA file
valid_genome_file = open('genome.txt', 'r')
valid_genome = valid_genome_file.readline()
assert(valid_genome == main.genome)

# Test reverse read
read = "ACGTTTANAGAC"
reversed_read = main.reverse_read(read)
assert(reversed_read == "CAGANATTTGCA")

# Test complement read
complemented_read = main.complement_read(read)
assert(complemented_read == "TGCAAATNTCTG")