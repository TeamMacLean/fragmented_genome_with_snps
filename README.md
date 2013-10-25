Fragmented Genome with SNPs
===========================

This is a model of a genome/chromosome/sequence from an individual created by the out-crossing of a mutant individual with non-mutants.
There is a high SNP density near the mutant position due to linkage.
The SNP density across this sequence follows a normal distribution, with the mutant position being located at 100,000 - the mid-point of this 200Kb sequence.
The sequence is fragmented, after being created, into fragments of 50-250b length. These model varying contig sizes.
A contig is a set of overlapping DNA segments that together represent a consensus region of DNA. My fragments are contigs, if they were real data, they would be formed by assembly of smaller overlapping reads.

### Making a fragmented genome with SNP's in Ruby:

1. [fragments_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/fragments_w_snps.rb) - This generates the fragments described above. JSON files are created: one contains the fragments, another the SNP positions for each fragment.

2. ^ a text file of the snp positions is also created. Running [normality_test.txt](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/normality_test.txt) in R confirms that the snp positions are normally distributed:
 - source("~/fragmented_genome_with_snps/normality_test.txt")
 - shapiro.test(snp_pos)

3. [json->fasta.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/json-%3Efasta.rb) - uses the information in the JSON files to construct a rudimentary fasta format file and VCF file.


### Designing an algorithm that will determine the position of a phenotype altering mutant, based on SNP density of the fragments:

1. [density_method_testa.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/density_method_testa.rb) - uses information from fasta and VCF files to work out the SNP density of each fragment in the fasta file, measured as SNPS per Kb.

2. Getting the fragments back into the correct order based on SNP density:
 - [rearrangement_methods.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/rearrangement_methods.rb) contains rearrangement methods
 - see rearrangement_methods.md for details of methods (and use of [munkres_test.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/munkres_test.rb))

3. rearrangement_methods.rb also creates text files that can be read by the R scripts ([random_vector_Rcode](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/random_vector_Rcode.txt) and [scatter_vectors](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/scatter_vectors.txt)) to create the graphs in [rearrangement_methods.md](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md)