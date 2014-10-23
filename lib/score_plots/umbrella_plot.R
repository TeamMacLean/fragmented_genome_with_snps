av_sd <- function(df, metric){
  names(df)[names(df)==metric] <- "metric_scores"
  df <- transform(df, Average = ave(metric_scores, gen, replicates, param_types))
  sd = aggregate(df$metric_scores, by=list(gen=df$gen, replicates=df$replicates, param_types=df$param_types), sd)
  df <- merge(df,sd)
  return(df)
}

uplot <- function(df, title, y_axis, metric, correct_fitness){
	df$param_types = factor(df$param_types, levels=c('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12'))
	df$replicates<- as.numeric(sub("p_run","",df[,2]))
	library(ggplot2)
	p <- ggplot(df, aes(colour = factor(replicates), y = metric_scores, x = gen)) +
	    xlab("Generations of genetic algorithm - each generation is a new population of contig permutations") +
	    ylab(y_axis) +
	    ggtitle(title) +
		scale_x_continuous() +
		geom_line(aes(y = Average), size=0.5) +
		geom_line(aes(y = Average-x), size=0.1, linetype="solid") +
		geom_line(aes(y = Average+x), size=0.1, linetype="solid") 
	if (correct_fitness != 'no_correct'){
		p <- p + geom_line(aes(y = correct_fitness, colour='C'), size=1.0, linetype="solid")
	}
	p <- p + facet_wrap(~param_types, ncol=2) +
	    theme_bw() +
	    theme(title = element_text(size = rel(0.5))) +
		guides(col = guide_legend(keywidth = 0.25, keyheight = 0.25, ncol = 4, byrow = TRUE, title.theme = element_text(size=8, angle = 0),
                              title = "Replicates"))
		if(metric!='Fitness'){
			p <- p + scale_y_continuous(limits=c(0, 1))
		}
	return(p)
}

av_sd2 <- function(df, metric){
  names(df)[names(df)==metric] <- "metric_scores"
  df <- transform(df, Average = ave(metric_scores, population))
  sd = aggregate(df$metric_scores, by=list(population=df$population), sd)
  df <- merge(df,sd)
  return(df)
}

metric_test_plot <- function(df, title, x_axis, y_axis, metric){
	df <- av_sd2(df, metric)
	library(ggplot2)
	p <- ggplot(df, aes(y = metric_scores, x = population)) +
	    xlab(x_axis) +
	    ylab(y_axis) +
	    ggtitle(title) +
	    # scale_y_continuous(limits=c(0, 1.2)) +
	    scale_y_continuous() +
		scale_x_continuous() +
		geom_line(aes(y = Average), size=0.5) + # Average of data
    	geom_line(aes(y = mean(df$shuffled[!is.na(df$shuffled)])), size=0.5) + # Average of random permutations
		geom_ribbon(aes(y = Average, ymin = (Average-x), ymax = (Average+x), fill="green"), alpha=0.5) + # SD of data
		geom_ribbon(aes(y = mean(df$shuffled[!is.na(df$shuffled)]), ymin = (mean(df$shuffled[!is.na(df$shuffled)]) - sd(df$shuffled[!is.na(df$shuffled)])),
                    ymax = (mean(df$shuffled[!is.na(df$shuffled)]) + sd(df$shuffled[!is.na(df$shuffled)])),
                    fill="red"), alpha=0.5) + # SD of random permutations
	    theme_bw() +
    guides(fill=FALSE) +
	  theme(title = element_text(size = rel(1.0)))
  return(p)
}