install.packages("UpSetR")
library("UpSetR")

Instance_count <- read.csv(file = "final_summary.txt", sep = "\t")

###Prepare the data file
names(Instance_count) <- sub("X", "", names(Instance_count))
row.names(Instance_count) <- Instance_count[,1]
Instance_count <- Instance_count[,-1]

###This converts all values to 0 or 1 (presence, absence)
Instance_count[] <- lapply(Instance_count, function(x) ifelse(x>0,1,x))

###Create plot
png(filename= "Per_species_identification_plot.png",width=min(3600,1800+800*ncol(plot)/10),height=1800,res=300)

upset(Instance_count, nsets = 20, order.by = c("freq", "degree"), main.bar.color="black", point.size = 3.5, line.size = 1, 
      mainbar.y.label = "Species intersection", sets.x.label = "Recovered sequence homologs")

dev.off()
