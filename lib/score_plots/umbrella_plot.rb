#encoding: utf-8
class UPlot
	require 'rinruby'
	require_relative '../write_it'
	require 'pdist'
	require_relative '../reform_ratio'

	# Input: The dataset to plot
	# Output: .csv file that can be loaded into R as a data frame, with the fitness scores and metric scores for each permutation
	def self.data_save(dataset)
		fasta = ReformRatio.fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")
		original_order = ReformRatio.fasta_id_n_lengths(fasta)[0]
		x = 1
		WriteIt.add_to("arabidopsis_datasets/#{dataset}/data.csv", 'gen,replicates,param_types,Fitness,Deviation,Square,Hamming,R,LCS,KT')
		Dir.entries("arabidopsis_datasets/#{dataset}").each do |run| 
			run_num = run.dup
			run_num.slice!('p_run')
			case run_num.to_i
				when 1..10 then param_type = 'p1' # Can be left the same for old and new experiment (old and new fitness methods)
				when 11..20 then param_type = 'p2'
				when 21..30 then param_type = 'p3'
				when 31..40 then param_type = 'p4'
				when 41..50 then param_type = 'p5'
				when 51..60 then param_type = 'p6'
				when 61..70 then param_type = 'p7'
				when 71..80 then param_type = 'p8'
				when 81..90 then param_type = 'p9'
				when 91..100 then param_type = 'p10'
				when 100..110 then param_type = 'p11'
				when 110..120 then param_type = 'p12'
			end
			if run.include?('p_run') # any other dir is not a replicate
				Dir.entries("arabidopsis_datasets/#{dataset}/#{run}").each do |gen| 
					if gen.include?('Gen')
						unless gen.include?('lists') || gen.include?('K')
							Dir.entries("arabidopsis_datasets/#{dataset}/#{run}/#{gen}").each do |ptxt|
								unless ptxt.include?('best') || ptxt.include?('table') # excluding table_data.txt and best_permutation.txt
									if ptxt.include?('.txt')
										perm = []
										IO.foreach("arabidopsis_datasets/#{dataset}/#{run}/#{gen}/#{ptxt}") { |line| perm << line.gsub(/\n/,'') }
										WriteIt.add_to("arabidopsis_datasets/#{dataset}/data.csv", 
											"#{gen[3..-1]},#{run},#{param_type},#{perm[0].to_f},#{PDist.deviation(original_order, perm[1..-1])},#{PDist.square(original_order, perm[1..-1])},#{PDist.hamming(original_order, perm[1..-1])},#{PDist.rdist(original_order, perm[1..-1])},#{PDist.lcs(original_order, perm[1..-1])},#{PDist.kendalls_tau(original_order, perm[1..-1])}")
										puts "permutation#{x}"
										x+=1
									end
								end
							end
						end
					end
				end
			end
		end
	end

	# Makes plot from arrays of generations (on for each data point), metric scores, and the run/replicate the data is from. Runs from different parameter replicates are faceted
	def self.uplot(dataset, filename, metric, y_axis, title, datafile, correct_fitness)
		myr = RinRuby.new(echo = false)
		myr.dataset = dataset
		myr.filename = filename
		myr.metric = metric
		myr.title = title
		myr.y_axis = y_axis
		myr.datafile = datafile
		myr.correct_fitness = correct_fitness
		myr.eval "source('~/fragmented_genome_with_snps/lib/score_plots/umbrella_plot.R')"
		myr.eval "df <- read.csv(file.path(paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset, sep=''), datafile))"
		myr.eval "df <- df <- av_sd(df, metric)"
		myr.eval "p <- uplot(df, title, y_axis, metric, correct_fitness)"
		myr.eval "ggsave(p, file = paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset,'/', filename,'.png', sep = ''))"
		myr.quit
		puts 'made a plot'
	end

	# Input: The dataset to plot
	# Output: .csv file that can be loaded into R as a data frame, with the fitness scores for each permutation
	def self.fits_save(dataset)
		fasta = ReformRatio.fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")
		original_order = ReformRatio.fasta_id_n_lengths(fasta)[0]
		x = 1
		WriteIt.add_to("arabidopsis_datasets/#{dataset}/data_fits.csv", 'gen,replicates,param_types,Fitness,Deviation,Square,Hamming,R,LCS,KT')
		Dir.entries("arabidopsis_datasets/#{dataset}").each do |run| 
			run_num = run.dup
			run_num.slice!('p_run')
			case run_num.to_i
				when 1..10 then param_type = 'p1' # Can be left the same for old and new experiment (old and new fitness methods)
				when 11..20 then param_type = 'p2'
				when 21..30 then param_type = 'p3'
				when 31..40 then param_type = 'p4'
				when 41..50 then param_type = 'p5'
				when 51..60 then param_type = 'p6'
				when 61..70 then param_type = 'p7'
				when 71..80 then param_type = 'p8'
				when 81..90 then param_type = 'p9'
				when 91..100 then param_type = 'p10'
				when 100..110 then param_type = 'p11'
				when 110..120 then param_type = 'p12'
			end
			if run.include?('p_run') # any other dir is not a replicate
				Dir.entries("arabidopsis_datasets/#{dataset}/#{run}").each do |gen| 
					if gen.include?('Gen')
						unless gen.include?('lists') || gen.include?('K')
							Dir.entries("arabidopsis_datasets/#{dataset}/#{run}/#{gen}").each do |ptxt|
								unless ptxt.include?('best') || ptxt.include?('table') # excluding table_data.txt and best_permutation.txt
									if ptxt.include?('.txt')
										perm = []
										IO.foreach("arabidopsis_datasets/#{dataset}/#{run}/#{gen}/#{ptxt}") { |line| perm << line.gsub(/\n/,'') }
										WriteIt.add_to("arabidopsis_datasets/#{dataset}/data_fits.csv", 
											"#{gen[3..-1]},#{run},#{param_type},#{perm[0].to_f}")
										puts "permutation#{x}"
										x+=1
									end
								end
							end
						end
					end
				end
			end
		end
	end
end

