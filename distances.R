#Analyse the genetic distance of the maverick compared to the surrounding 
# non-coding regions

library(ape) #Load ape
library(ggplot2)

#Load fasta MSAs of the corresponding regions using ape

size7 = 123980-110963+1 #maverick size in chr7 alignment
size8 = 142050-141408+1 #maverick size in chr8 alignment

start7 = 110963 #start position of the chr 7 maverick in the alignment 
start8 = 141408 #start position of the chr8 maverick in the alignment

align7 = 404032 #alignment size in base pairs for chr7
align8 = 315835 #alignment size in base pairs for chr8

chr7_1 = c(1,start7-size7-1) #upstream sampling interval
chr7_2 = c(start7+size7+1,align7-size7) #downstream sampling interval

chr8_1 = c(1,start8-size8-1) #upstream sampling interval
chr8_2 = c(start8+size8+1,align8-size8) #downstream sampling interval

#Define function to obatin coordinates of samplig intervals

sample_alignment <- function(coordinates,n,size) {
  
  samples <- sample(coordinates[1]:coordinates[2],n)
  
  df <- data.frame(Pos1=integer(), Pos2=integer())
  
  for (i in samples) {
    
    df <- rbind(df,data.frame(Pos1=i,Pos2=i+size))
    
  } 
  
  return(df)       
  
}

# Read MSAs from file

chr7.dna <- read.dna("chr7-region-revcomp.fasta", format = "fasta")
chr8.dna <- read.dna("chr8-region-revcomp.fasta", format = "fasta")

#Define distance function (observed genetic distances as percentages)

genetic_distances <- function(dataframe,alignment) { #takes in data frame of alignment coordinates
  
  distances <- c()
  
  for (i in 1:nrow(dataframe)){
    
    start = dataframe[i,1]
    end = dataframe[i,2]
    
    subset = alignment[1:dim(alignment)[1],1:dim(alignment)[2]]
    
    distances <- c(distances,as.vector(dist.gene(subset,method="percentage")))
    
    cat('Row', i, "done!",'\n')
    
  }
  
  return(distances)
  
}

#Calculate genetic distances outside mavericks

chr7_samples <- rbind(sample_alignment(chr7_1,50,size7),sample_alignment(chr7_2,50,size7))
chr7_distances <- genetic_distances(chr7_samples,chr7.dna)
write.csv(chr7_distances,file="chr7_distances.csv")

chr8_samples <- rbind(sample_alignment(chr8_1,50,size8),sample_alignment(chr8_2,50,size8))
chr8_distances <- genetic_distances(chr8_samples,chr8.dna)
write.csv(chr8_distances,file="chr8_distances.csv")

#Calculate genetic distances of mavericks

chr7_maverick <- data.frame(Pos1=110963,Pos2=123980)
chr7_maverick_distance <- genetic_distances(chr7_maverick,chr7.dna)

chr8_maverick <- data.frame(Pos1=141408,Pos2=142050)
chr8_maverick_distance <- genetic_distances(chr8_maverick,chr8.dna)

#Plot histograms and overlay distributions

chr7_noncoding <- data.frame(dna=as.factor(rep("Non-coding",length(chr7_distances))),distances=chr7_distances)
chr7_maverick <- data.frame(dna=as.factor(rep("Maverick",length(chr7_maverick_distance))),distances=chr7_maverick_distance)
chr7 <- rbind(chr7_noncoding,chr7_maverick)

chr8_noncoding <- data.frame(dna=as.factor(rep("Non-coding",length(chr8_distances))),distances=chr8_distances)
chr8_maverick <- data.frame(dna=as.factor(rep("Maverick",length(chr8_maverick_distance))),distances=chr8_maverick_distance)
chr8 <- rbind(chr8_noncoding,chr8_maverick)

ggplot(chr7, aes(x=distances,fill=dna)) +
  geom_density(alpha=0.4) +
  xlab("Observed genetic distance (proportion)") +
  ylab("Density") +
  scale_fill_discrete(name = "") +
  theme_bw() +
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        legend.text=element_text(size=12))

ggplot(chr8, aes(x=distances,fill=dna)) +
  geom_density(alpha=0.4) +
  xlab("Observed genetic distance (proportion)") +
  ylab("Density") +
  scale_fill_discrete(name = "") +
  theme_bw() +
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        legend.text=element_text(size=12))

# Kolmogorov-Smirnov test for the equality of the distributions

ks.test(chr7_maverick$distances, chr7_noncoding$distances, alternative = c("less"))

ks.test(chr8_maverick$distances, chr8_noncoding$distances, alternative = c("less"))
