#encoding: utf-8
class ExamplePerms

	# Input 0: Array of fasta format fragments
	# Input 1: Permutation array of frag ids
	# Output: Array fasta frags re-ordered by the permutation ids
	def self.fasta_p_id(fasta, perm_ids)
		fasta_perm = []
		perm_ids.each do |id|
			fasta.each do |frag|
				if id == frag.entry_id
					fasta_perm << frag
				end
			end
		end
		return fasta_perm
	end

	# Input 0: Number of generations to plot
	# Input 1: The generation to start with
	# Input 2: The number of generations to increment by
	# Input 3: Dataset
	# Input 4: Run
	# Output: Array of populations, each population is from the next generation and has a number of permutations (arrays of fasta frag ids), with the fitness score at element 0 of each
	def self.get_perms(gen, start, inc, dataset, run) # number of generations to choose from
		all_perms = []
		n = start
		gen.times do
			pop = []
			Dir.entries("arabidopsis_datasets/#{dataset}/#{run}/Gen#{n}").each do |ptxt|
				unless ptxt.include?('best') || ptxt.include?('table') # excluding table_data.txt and best_permutation.txt
					if ptxt.include?('.txt')
						perm = []
						IO.foreach("arabidopsis_datasets/#{dataset}/#{run}/Gen#{n}/#{ptxt}") { |line| perm << line.gsub(/\n/,'') }
						pop << perm
					end
				end
			end
			all_perms << pop
			n+=inc
		end
		return all_perms
	end
end

