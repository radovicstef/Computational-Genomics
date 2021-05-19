import os

import main


print("Testing started...")


# Run simulate method with simple sequence sample
main.simulate("genomeSample.fa", 60, 10, 4, 7)


# Tests


# Test if sample sequence is correctly extracted from FASTA file as one string
# Valid genome sequence should be in genomeSample.txt in one line
valid_genome_file = open('genomeSample.txt', 'r')
valid_genome = {}
genome_name = ''
lines = valid_genome_file.read().split()
for line in lines:
    if line[0] == '>':
        genome_name = line[1::]
        valid_genome[genome_name] = ''
    else:
        valid_genome[genome_name] = line
for sequence_name, sequence in main.genome.items():
    assert sequence == valid_genome[sequence_name]


# Sample reads
reads = ["ACCTAG", "ACACGGGTACC", "ACCCAGTTTTACGT", "ACGACG"]


# Test reverse read
reverse_reads = ["GATCCA", "CCATGGGCACA", "TGCATTTTGACCCA", "GCAGCA"]
reversed_read_generated = [main.reverse_read(read) for read in reads]
for i in range(len(reads)):
    assert(reverse_reads[i] == reversed_read_generated[i])


# Test complement read
complemented_reads = ["TGGATC", "TGTGCCCATGG", "TGGGTCAAAATGCA", "TGCTGC"]
complemented_reads_generated = [main.complement_read(read) for read in reads]
for i in range(len(reads)):
    assert(complemented_reads[i] == complemented_reads_generated[i])

# Test get quality
# The function should return read_length number of single nucleotide qualities based on average value
avg_quality = 15
sigma = 5
read_length = 5
qualities = main.get_quality(avg_quality, sigma, read_length)
for quality in qualities:
    assert (quality >= chr(avg_quality - sigma) or quality <= chr(avg_quality + sigma))


# Test getting the mutated single nucleotide
assert main.get_mutated_single_nucleotide("A") != "A"
assert main.get_mutated_single_nucleotide("C") != "C"
assert main.get_mutated_single_nucleotide("G") != "G"
assert main.get_mutated_single_nucleotide("T") != "T"
assert main.get_mutated_single_nucleotide("N") != "N"


# Test the SNV mutations
# Number of mutations can be less or equal to number of expected SNV mutations
# as it is possible to mutate nucleotide more than once
prob_snv = 0.2
number_of_snv_mutations_expected = round(len(main.genome[sequence_name]) * prob_snv)
main.add_mutations(prob_snv)
for sequence_name, sequence in main.genome.items():
    number_of_snv_mutations = 0
    for i in range(len(valid_genome[sequence_name])):
        if sequence[i] != valid_genome[sequence_name][i]:
            number_of_snv_mutations += 1
    assert(number_of_snv_mutations <= number_of_snv_mutations_expected)


# Test the INS mutations
prob_ins = 0.2
number_of_ins_mutations_expected = {}
for sequence_name, sequence in main.genome.items():
    number_of_ins_mutations_expected[sequence_name] = round(len(main.genome[sequence_name]) * prob_ins)
main.add_mutations(0, prob_ins)
for sequence_name, sequence in main.genome.items():
    assert(number_of_ins_mutations_expected[sequence_name] == len(main.genome[sequence_name]) - len(valid_genome[sequence_name]))

# Test the DEL mutations
prob_del = 0.2
print("Genome length " + str(len(main.genome[sequence_name])))
genome_current_length = len(main.genome[sequence_name])
number_of_del_mutations_expected = {}
for sequence_name, sequence in main.genome.items():
    number_of_del_mutations_expected[sequence_name] = round(len(main.genome[sequence_name]) * prob_del)
main.add_mutations(0, 0, prob_del)
for sequence_name, sequence in main.genome.items():
    assert(number_of_del_mutations_expected[sequence_name] == genome_current_length - len(main.genome[sequence_name]))


# Test input parameters
assert (main.check_input_parameters(60, 10, 12, 24) == 0)


print("Finished!")

