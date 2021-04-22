import main

# Run main method
main.simulate("genomeSample.fa", 100, 4, 7, 10)

# Tests

# Test if sample genome is correctly extracted from FASTA file
valid_genome_file = open('genome.txt', 'r')
valid_genome = [line.rstrip('\n') for line in valid_genome_file.readlines()]
i = 0
for sequence_name, sequence in main.genome.items():
    assert sequence == valid_genome[i]
    i += 1

# Test reverse read
read = "ACGTTTANAGAC"
reversed_read = main.reverse_read(read)
assert(reversed_read == "CAGANATTTGCA")

# Test complement read
complemented_read = main.complement_read(read)
assert(complemented_read == "TGCAAATNTCTG")

# Test get quality
avg_quality = 15
sigma = 5
read_length = 5
qualities = main.get_quality(avg_quality, sigma, read_length)
print(qualities)
for quality in qualities:
    assert (quality >= chr(avg_quality - sigma) or quality <= chr(avg_quality + sigma))


# Test getting the mutated single nucleotide
assert main.get_mutated_single_nucleotide("A") != "A"
assert main.get_mutated_single_nucleotide("C") != "C"
assert main.get_mutated_single_nucleotide("G") != "G"
assert main.get_mutated_single_nucleotide("T") != "T"
assert main.get_mutated_single_nucleotide("N") != "N"


# Test the mutations
main.add_mutations(0.2)
print("MUTATED GENOME: ")
for sequence_name in main.genome:
    print(main.genome[sequence_name])
