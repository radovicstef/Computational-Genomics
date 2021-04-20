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

# Test get quality
avg_quality = 15
sigma = 5
read_length = 5
qualities = main.get_quality(avg_quality, sigma, read_length)
print(qualities)
for quality in qualities:
    assert (quality >= chr(avg_quality - sigma) or quality <= chr(avg_quality + sigma))
