#calculate differences in introns
setwd("/Users/rishidek/Dropbox/RishiMAC/Danaus/Genome_annotation/introns/auto_intron_exon/t1_only/")

dchr2_2_EXONS <- read.csv("Dchry2_2.exons_length_output.txt", head = FALSE)
dchr2_2_INTRONS <- read.csv("Dchry2_2.introns_length_output.txt", head = FALSE)

dplex_V4_B_EXONS <- read.csv("dplex_v4_B.exons_length_output.txt", head = FALSE)
dplex_V4_B_INTRONS <- read.csv("dplex_v4_B.introns_length_output.txt", head = FALSE)

dplex_mex_B_EXONS <- read.csv("dplex_mex_B.exons_length_output.txt", head = FALSE)
dplex_mex_B_INTRONS <- read.csv("dplex_mex_B.introns_length_output.txt", head = FALSE)

dchr2_2_EXONS <- subset(dchr2_2_EXONS, dchr2_2_EXONS$V1 > 0)
dchr2_2_INTRONS <- subset(dchr2_2_INTRONS, dchr2_2_INTRONS$V1 > 0)

dplex_V4_B_EXONS <- subset(dplex_V4_B_EXONS, dplex_V4_B_EXONS$V1 > 0)
dplex_V4_B_INTRONS <- subset(dplex_V4_B_INTRONS, dplex_V4_B_INTRONS$V1 > 0)

dplex_mex_B_EXONS <- subset(dplex_mex_B_EXONS, dplex_mex_B_EXONS$V1 > 0)
dplex_mex_B_INTRONS  <- subset(dplex_mex_B_INTRONS, dplex_mex_B_INTRONS$V1 > 0)


mean(dchr2_2_EXONS$V1)
mean(dchr2_2_INTRONS$V1)

mean(dplex_V4_B_EXONS$V1)
mean(dplex_V4_B_INTRONS$V1)

mean(dplex_mex_B_EXONS$V1)
mean(dplex_mex_B_INTRONS$V1)

par(mfrow=c(1,2))
boxplot(dchr2_2_EXONS$V1, dplex_mex_B_EXONS$V1, dplex_V4_B_EXONS$V1, xaxt = 'n', main = "Exons", ylab = "bp", outline = F, col = c("goldenrod", "tomato", "turquoise4"))
axis(side = 1, at = c(1:3), labels = c(paste("Dchr2.2 - mean=", as.integer(round(mean(dchr2_2_EXONS$V1))), sep = ""),
                                       paste("MEX_DaPlex - R - mean=", as.integer(round(mean(dplex_mex_B_EXONS$V1))), sep = ""), 
                                       paste("Dplex_V4 - R - mean=", as.integer(round(mean(dplex_V4_B_EXONS$V1))), sep = "")), cex.axis = 0.5)

boxplot(dchr2_2_INTRONS$V1, dplex_mex_B_INTRONS$V1,  dplex_V4_B_INTRONS$V1, xaxt = 'n', main = "Introns", ylab = "bp", outline = F, col = c("goldenrod", "tomato", "turquoise4"))
axis(side = 1, at = c(1:3), labels = c(paste("Dchr2.2 - mean=", as.integer(round(mean(dchr2_2_INTRONS$V1))), sep = ""), 
                                       paste("MEX_DaPlex - R - mean=", as.integer(round(mean(dplex_mex_B_INTRONS$V1))), sep = ""),
                                       paste("Dplex_V4 - R - mean=", as.integer(round(mean(dplex_V4_B_INTRONS$V1))), sep = "")), cex.axis = 0.5)

?wilcox.test
wilcox.test(dchr2_2_INTRONS$V1, dplex_mex_B_INTRONS$V1, c("greater"))
wilcox.test(dchr2_2_INTRONS$V1, dplex_V4_B_INTRONS$V1, c("greater"))
wilcox.test(dplex_mex_B_INTRONS$V1, dplex_V4_B_INTRONS$V1, c("less"))

wilcox.test(dchr2_2_EXONS$V1, dplex_mex_B_EXONS$V1, c("greater"))
wilcox.test(dchr2_2_EXONS$V1, dplex_V4_B_EXONS$V1, c("greater"))
wilcox.test(dplex_mex_B_INTRONS$V1, dplex_V4_B_INTRONS$V1, c("less"))


val <- c(mean(dchr2_2_EXONS$V1), mean(dchr2_2_INTRONS$V1), "red")
val3 <- c(mean(dplex_V4_B_EXONS$V1), mean(dplex_V4_B_INTRONS$V1), "lightblue")
val5 <- c(mean(dplex_mex_B_EXONS$V1), mean(dplex_mex_B_INTRONS$V1), "lightgreen")

all_vals <- as.data.frame(rbind(val, val3, val5))
plot(as.integer(as.character(all_vals$V1)) ~ as.integer(as.character(all_vals$V2)), col = as.character(all_vals$V3), xlim = c(500,1100), ylab = "exon length - bp", xlab = "intron length - bp")


total_introns_dchr2_2 <- sum(dchr2_2_INTRONS$V1)
total_introns_dchr2_2
#80067597
total_introns_dplexmex_B <- sum(dplex_mex_B_INTRONS$V1)
total_introns_dplexmex_B
#60058038
total_introns_dplex_v4_B <- sum(dplex_V4_B_INTRONS$V1)
total_introns_dplex_v4_B
#64553725

list_introns <- c(total_introns_dchr2_2, total_introns_dplexmex_B, total_introns_dplex_v4_B)
numb <- c(1:3)
plot(numb, list_introns, ylim = c(0, max(list_introns)))

total_introns_dchr2_2-total_introns_dplexmex_B
#20009559
total_introns_dchr2_2-total_introns_dplex_v4_B
#15513872

