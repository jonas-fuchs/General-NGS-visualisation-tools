##library import
library(ggplot2)
library(dplyr)
library(haven)
library(scales)
library(tools)
library(fs)
library(data.table)
library(phylotools)
library(base)
library(stringr)

##settings
min_coverage <- 20

##read in reference fasta and determine name and length
reference_fasta <- read.fasta(file = "reference.fasta")
sequences <- data.table()
sequences$seq_name <- str_extract(reference_fasta$seq.name, pattern = "[:graph:]{1,}")
sequences$seq_length <- nchar(reference_fasta$seq.text, type = "chars")
sequences <- sequences[order(sequences$seq_length, decreasing = T),]


## define limits for the plots as 3' and 5' end might be not sufficiently covered
limit_for_plot <- data.table()
limit_for_plot$V1 <- sequences$seq_name
limit_for_plot$V2 <- 0
limit_for_plot$V3 <- 1
limit_for_plot_2 <- data.table()
limit_for_plot_2$V1 <- sequences$seq_name
limit_for_plot_2$V2 <- sequences$seq_length+1
limit_for_plot_2$V3 <- 1

limit_for_plot <- rbind(limit_for_plot,limit_for_plot_2)


##get input files
if(file.exists("Tabular_coverage_per_sequence.tsv")){
        file.remove("Tabular_coverage_per_sequence.tsv")
}

if(file.exists("Tabular_coverage_mean.tsv")){
        file.remove("Tabular_coverage_mean.tsv")
}


tsv_files <- list.files(path = ".", pattern = "*.tsv",recursive = TRUE ,full.names = TRUE)

tabular_overview <- data.table()
tabular_overview_mean <- data.table()

pdf(paste0("Coverage", ".pdf"), width = 15, height = 10)

for (i in 1:length(tsv_files)) {
        
        tsv_files_temp <- read.table(tsv_files[i], sep = "\t")
        tsv_files_temp <- tsv_files_temp[order(match(tsv_files_temp$V1, sequences$seq_name)),]
        
        Sample_name <- path_file(tsv_files[i])
        Sample_name = file_path_sans_ext(Sample_name)

        tabular_overview_temp <- data.table()
        tabular_overview_mean_temp <- data.table()
        
        for (j in 1:length(unique(tsv_files_temp$V1))) {
                tsv_split_temp <- tsv_files_temp[tsv_files_temp$V1 == unique(tsv_files_temp$V1)[j],]
                
                ## get non covered regions
                tabular_overview_temp_2 <- data.table()
                
                seq_length_temp <- as.numeric(sequences$seq_length[sequences$seq_name == unique(tsv_files_temp$V1)[j]])
                non_covered <- as.numeric(length(which(tsv_split_temp$V3 < min_coverage))) + seq_length_temp - as.numeric(length(tsv_split_temp$V3))
                
                tabular_overview_temp_2$Sample_name <- Sample_name
                seq_name <- sequences$seq_name[sequences$seq_name == unique(tsv_files_temp$V1)[j]]
                tabular_overview_temp_2$Seq_name <- seq_name
                tabular_overview_temp_2$Percent_coverage <- format(round(100 - (non_covered/seq_length_temp*100)), nsmall = 2)
                
                #mean coverage
                mean_coverage <- sum(tsv_split_temp$V3)/seq_length_temp
                mean_coverage <- format(round(mean_coverage,0), nsmall = 0)
                tabular_overview_temp_2$Mean_coverage <- mean_coverage
                
                tabular_overview_temp <- rbind(tabular_overview_temp, tabular_overview_temp_2)
                
        }
        
        ## get coverage and mean coverage over the whole genome
        percent_covered <- sum(as.numeric(tabular_overview_temp$Percent_coverage))/nrow(tabular_overview_temp)
        percent_covered <- format(round(percent_covered,0), nsmall = 2)
        mean_coverage_total <- sum(as.numeric(tabular_overview_temp$Mean_coverage))/nrow(tabular_overview_temp)
        mean_coverage_total <- format(round(mean_coverage_total,0), nsmall = 0)
        
        tabular_overview_mean_temp$Sample_name <- Sample_name
        tabular_overview_mean_temp$Percent_coverage <- percent_covered
        tabular_overview_mean_temp$Mean_coverage <- mean_coverage_total
        
        tabular_overview_mean <- rbind(tabular_overview_mean, tabular_overview_mean_temp)
        
        ## get values for mean and min for plot
        dMean <- tsv_files_temp %>%
                group_by(V1) %>%
                summarise(MN = mean(V3))
        dMin <- tsv_files_temp %>%
                group_by(V1) %>%
                summarise(MN = min(V3))
        
        
        ## reformat tsv data
        tsv_files_temp <- as.data.table(tsv_files_temp)
        tsv_files_temp <- rbind(tsv_files_temp, limit_for_plot) ## bind plot limits to table
        tsv_files_temp$V1 <- factor(tsv_files_temp$V1, levels = unique(tsv_files_temp$V1)) ## order plots based on sequence length
        
        #define upper limit for y-axis
        upper_limit <- 10^ceiling(log10(max(tsv_files_temp$V3)))
        
        ## plot data in facet plot
        plot <- ggplot(tsv_files_temp, aes(x=V2,
                                           y=V3)) +
                geom_area(color="gray55",
                          fill= "lightsteelblue",
                          alpha=0.6) +
                scale_x_continuous(name = "nucleotide position",
                                   expand = expansion(mult = c(0,0))) +
                scale_y_continuous(name = "Coverage",
                                   trans = "log10", limits = c(1,upper_limit),
                                   expand = expansion(mult = c(0,0)),
                                   breaks = trans_breaks("log10", function(x) 10^x), 
                                   labels = trans_format("log10", math_format(10^.x))) +
                annotation_logticks(sides = "l",
                                    size = 0.2,
                                    outside = TRUE) +
                coord_cartesian(clip="off") +
                facet_wrap(~ V1,
                           scales = "free") +
                geom_hline(yintercept=min_coverage,
                           linetype="dashed",
                           color = "red",
                           size = 0.2) +
                geom_text(aes(x = 10, y = min_coverage),
                          label = "20x",
                          size = 2,
                          vjust = 1,
                          hjust = 0,
                          color = "red") +
                geom_hline(data = dMean,
                           aes(yintercept = MN),
                           linetype = "dotted",
                           color = "grey40",
                           size = 0.2)+
                geom_text(data = dMean, aes(x = 10, y = MN),
                          label = "average",
                          size = 2,
                          vjust = 1,
                          hjust = 0,
                          color = "grey 40") +
                geom_hline(data = dMin,
                           aes(yintercept = MN),
                           linetype = "dashed",
                           size = 0.2) +
                geom_text(data = dMin,
                          aes(x = 10, y = MN),
                          label = "min",
                          size = 2,
                          vjust = 1,
                          hjust = 0) +
                labs(title = paste0(Sample_name),
                     subtitle = paste0("% covered = ",percent_covered, ", ", "Overall mean coverage = ", mean_coverage_total, "x")) +
                theme(axis.line = element_line(color = "black", size = 0.2),
                      axis.ticks = element_line(colour = "black", size = 0.2),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      axis.text.x = element_text(face = "plain", color = "black", size = 9, hjust = 0.5),
                      axis.text.y = element_text(face = "plain", color = "black", size = 9, hjust = 0.5, margin = margin(r = 8)),
                      plot.title = element_text(face = "bold", color = "black", size = 12, hjust = (0.5)),
                      plot.margin = unit(c(1,1,1,1), "cm")
                 )
                
        print(plot)
        tabular_overview <- rbind(tabular_overview, tabular_overview_temp)
}

dev.off()

if(nrow(sequences) > 1){
        fwrite(tabular_overview, "Tabular_coverage_per_sequence.tsv", sep = "\t")  
}

fwrite(tabular_overview_mean, "Tabular_coverage_mean.tsv", sep = "\t")

