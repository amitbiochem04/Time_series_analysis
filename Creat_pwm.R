
#########################
# load/install libraries
#########################

# better thatn: source("http://bioconductor.org/biocLite.R") is nowadays:
if (!require(BiocManager)) {
  install.packages("BiocManager")
  require(BiocManager)
}

if (!require(MatrixRider)) {
  BiocManager::install("MatrixRider")
  require(MatrixRider)
}

if (!require(Biostrings)) {
  BiocManager::install("Biostrings")
  require(Biostrings)
}

##################################
# create a position weight matrix
##################################

dna2consensusMatrix <- function(dna_sequence) {
  dna = DNASTringSet(dna_equence)
  consensusMatrix(dna)
}

dna2pwm <- function(dna_sequence) {
  PWM(dna2consensusMatrix(dna_sequence))
}

pwm <- dna2pwm("ATTCCATTGATCACGA")



# pwm is your position weight matrix
# now create a poistion profile matrix to get a score 

############################################################
# Looking for binding potential for a single TF on a sequence
############################################################

d <- dna2consensusMatrix(dna_sequence)[1:4,]

tf_bindingsites <- function(search_seq, tf_seq, 
                            ID="LRE", 
                            name="FRQ",
                            matrixClass="Unknown",
                            strand="+",
                            tags=list(),
                            bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                           ...) {
  seq <- unlist(readDNAStringSet("selected_gene_seq.fa"))
  pfm <- PFMatrix(
                profileMatrix=d)
  getSeqOccupancy(sequence, pfm, 1)
}



pfm <- PFMatrix(ID="LRE", 
                name="FRQ",
                matrixClass="Unknown",
                strand="+",
                tags=list(),
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                profileMatrix=d, ...)

# load the fasta file, meaning that gene of interest here we have FRQ

sequence <- readDNAStringSet("selected_gene_seq.fa")
sequence <- unlist(sequence)

# the sequence should be in dnastring format 
# sequence <- DNAString("TACCTACGTCTATGACTCCATCGGCCTCTCTTTCCTCCCTTCGGCGCAACAGAGAGCGTTCATCAAAGGGCTGATCAATGATTTAAAG")

getSeqOccupancy(sequence, pfm, 1)

# library(MatrixRider) only work in new verison (R 3.4 version) so you can update 
# below is some code which will help you to restore all your R package even if you upgrade your R 

## Change accordingly
list_dir <- "/Users/amit/Desktop/Axel/"

## Get the list of installed packages
installed <- dir(.libPaths())

## Save the list for later use
save(installed, 
     file = file.path(list_dir, 
                      paste0(Sys.Date(), 
                             '-installed.Rdata')))

## Explore the list
head(installed)

## [1] "abind"   "acepack" "ada"     "AER"     "affy"    "affyio"

length(installed)

# now install new verison of R from CRAN,
# then follow the command once you install open R 
# then follow the command 

## Change accordingly
list_dir <- "/Users/amit/Desktop/Axel/"

## Find the corresponding Rdata files
previous <- dir(path = list_dir, pattern = 'installed.Rdata')

## Load the latest one
load(file.path(list_dir, previous[length(previous)]))

## Just checking it
head(installed)

## [1] "abind"   "acepack" "ada"     "AER"     "affy"    "affyio"
#Next, get the list of current R packages you have installed. Every new R installation comes with a few of them (the base packages). You don’t need to install those.

current <- dir(.libPaths())

# Finally, install the missing packages
# source('http://bioconductor.org/biocLite.R')
# biocLite(installed[!installed %in% current])
require(BiocManager)
BiocManager::install(installed[!installed %in% current])

# Check which packages are missing
current_post_installation <- dir(.libPaths())
installed[!installed %in% current_post_installation]

###install from github
if (!require(devtools)) {
  install.packages('devtools')
  require('devtools')
}
devtools::install_github("jalvesaq/colorout")

