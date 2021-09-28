#repeat landscapes of danaus genomes

repeats_df <- read.csv(file = "../../repeats/repeat_summaries.csv")

length_dchr <- 354020300
length_dplex4 <- 248703116
length_dplexM <- 245173502

Dchry <- as.data.frame(repeats_df$Dchr2_2[4:8])
Dple4 <- as.data.frame(repeats_df$Dplex_V4[4:8])
DpleM <- as.data.frame(repeats_df$Dplex_mex[4:8])

Dchry$percent <- (Dchry$`repeats_df$Dchr2_2[4:8]`/length_dchr)*100
Dple4$percent <- (Dple4$`repeats_df$Dplex_V4[4:8]`/length_dplex4)*100
DpleM$percent <- (DpleM$`repeats_df$Dplex_mex[4:8]`/length_dplexM)*100

Dchry$numb <- 1:5
Dple4$numb <- 1:5
DpleM$numb <- 1:5

Dchry$spp <- "Dchry"
Dple4$spp <- "Dplex_zV4"
DpleM$spp <- "Dplex_mex"

fams <- c("masked", "retroelements", "DNA_transposons",  "rolling_circles", "unclassified")

Dchry$fams <- fams
Dple4$fams <- fams
DpleM$fams <- fams

a <- colnames(Dchry)[2:5]
b <- c("bp", a)

colnames(Dchry) <- b
colnames(Dple4) <- b
colnames(DpleM) <- b

full <- rbind(Dchry, Dple4, DpleM)

full_ordered <- full[order(full$spp),]
full_ordered <- full_ordered[order(full_ordered$numb),]

unique(full_ordered$fams)

masked <- subset(full_ordered, full_ordered$fams=="masked")
retroelements <- subset(full_ordered, full_ordered$fams=="retroelements")
DNA_transposons <- subset(full_ordered, full_ordered$fams=="DNA_transposons")
rolling_circles <- subset(full_ordered, full_ordered$fams=="rolling_circles")
unclassified <- subset(full_ordered, full_ordered$fams=="unclassified")

par(mfrow=c(1,5))
par(mar=c(4,2,4,0.2))
barplot(masked$percent, col=c("goldenrod","tomato", "turquoise4"), ylim = c(0,40), xlab = "")
title(xlab="% masked", line=1, cex.lab=1.2)
barplot(retroelements$percent, col=c("goldenrod","tomato", "turquoise4"), ylim = c(0,40), yaxt = 'n', xlab = "")
title(xlab="% Retroelements", line=1, cex.lab=1.2)
barplot(DNA_transposons$percent, col=c("goldenrod","tomato", "turquoise4"), ylim = c(0,40), yaxt = 'n', xlab = "")
title(xlab="% DNA transposons", line=1, cex.lab=1.2)
barplot(rolling_circles$percent, col=c("goldenrod","tomato", "turquoise4"), ylim = c(0,40), yaxt = 'n', xlab = "")
title(xlab="% Rolling-circles", line=1, cex.lab=1.2)
barplot(unclassified$percent, col=c("goldenrod","tomato", "turquoise4"), ylim = c(0,40), yaxt = 'n', xlab = "")
title(xlab="% Unclassified", line=1, cex.lab=1.2)
legend("topright", legend=c("Dchr2.2", "MEX_DaPlex", "Dplex_v4"),
       col=c("goldenrod", "tomato", "turquoise4"), pch = 15, cex=1)


