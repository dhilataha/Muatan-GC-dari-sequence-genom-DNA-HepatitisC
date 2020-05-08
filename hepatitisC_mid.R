#Fadhillah Putri Taha H071171301
#Gene Bank Data : Hepatitis C -- Hepatitis C virus full-length replicon pFGR-JFH1 RNA, complete sequence
#Dengan kode NCBI AB237837.1
#https://www.ncbi.nlm.nih.gov/nuccore/AB237837.1

#2a
library("seqinr")
hepatitisC <- read.fasta(file = "D:/dhil/Tugas/Bioinformatika/Bioinformatics/MID/hepatitisC.fasta")
hepatitisC_seq <- hepatitisC[[1]]

#2b
length(hepatitisC_seq)

#2c
table(hepatitisC_seq)

#2d
#Kadar konten GC
#Content GC = (Jumlah G = jumlah C)*100/panjang genom
((3127 + 3315)*100)/(2255 + 3315 + 3127 + 2414)
GC(hepatitisC_seq)

#2e
table(hepatitisC_seq)

#2f
count(hepatitisC_seq, 2)

#2g
first_data <- head(hepatitisC_seq,1000)
count(first_data, 2)

last_data <- tail(hepatitisC_seq, 1000)
count(last_data, 2)




# B_1 & B2
slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  hepatitisC_GCs <- numeric(n)
  for (i in 1:n) {
    hepatitisC <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    hepatitisC_GC <- GC(hepatitisC)
    hepatitisC_GCs[i] <- hepatitisC_GC
  }
  plot(starts,hepatitisC_GCs,type="b",xlab="Nucleotide start position",ylab="GC content")
  cat("\npuncak= ", max(hepatitisC_GCs))
  cat("\npalung= ", min(hepatitisC_GCs))
}

#a sliding window plot with a window size of 200 nucleotides
slidingwindowplot(200, hepatitisC_seq)
#a sliding window plot with a window size of 20000 nucleotides
slidingwindowplot(20000, hepatitisC_seq)
#a sliding window plot with a window size of 200000 nucleotides
slidingwindowplot(200000, hepatitisC_seq)

#a sliding window plot with a window size of 50 nucleotides
slidingwindowplot(50, hepatitisC_seq)
#a sliding window plot with a window size of 5000 nucleotides
slidingwindowplot(5000, hepatitisC_seq)
#a sliding window plot with a window size of 11000 nucleotides
slidingwindowplot(11000, hepatitisC_seq)

#a sliding window plot with a window size of 6000 nucleotides
slidingwindowplot(6000, hepatitisC_seq)

#B_3
AT <- function(inputseq)
{
  mytable <- count(inputseq, 1) # make a table with the count of As, Cs, Ts, and Gs
  mylength <- length(inputseq) # find the length of the whole sequence
  myAs <- mytable[[1]] # number of As in the sequence
  myTs <- mytable[[4]] # number of Ts in the sequence
  myAT <- (myAs + myTs)/mylength
  return(myAT)
}
AT(hepatitisC_seq)
GC(hepatitisC_seq)
0.4202142 + 0.5797858

#B_4
slidingwindowplotAT <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  hepatitisC_ATs <- numeric(n)
  for (i in 1:n) {
    hepatitisC <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    heppatitisC_AT <- AT(hepatitisC)
    hepatitisC_ATs[i] <- hepatitisC_AT
  }
  plot(starts,hepatitisC_ATs,type="b",xlab="Nucleotide start position",ylab="AT content")
  cat("\npuncak= ", max(hepatitisC_ATs))
  cat("\npalung= ", min(hepatitisC_Ats))
}
#a sliding window plot with a window size of 200 nucleotides
slidingwindowplot(200, hepatitisC_seq)

#B_5
count(hepatitisC_seq, 3)
sum(count(hepatitisC_seq,3))
199/11109

#B_5c
count(hepatitisC_seq,1)
2255/11111 #frekuensi A
3315/11111 #frekuensi C
3127/11111 #frekuensi G
0.202952 * 0.298353 * 0.2814328 #frekuensi GAC

#menghitung nilai rho
(199/11109)/(0.202952 * 0.298353 * 0.2814328)

#menghitung panjang rho untuk urutan 3 nukleotida
rho(hepatitisC_seq, 3)

#5d
help.search("rho")
base::getHook
seqinr::rho
