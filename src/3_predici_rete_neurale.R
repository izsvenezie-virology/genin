#dati iniziali
options <- commandArgs(trailingOnly = TRUE)
neural_net_file <- as.character(options[1]) #file contenente l'immagine della rete neurale prodotta
output <- as.character(options[2]) #file di output nel quale stampare le probabilità predette per la sequenza ignota

#carico le librerie necessari
library("nnet")

#prendo i dati da standard input
raw_data <- read.csv("stdin", header=FALSE, sep="\t", na.strings="", stringsAsFactors=FALSE, check.names=FALSE)

#carico la rete neurale pre-prodotta
load(neural_net_file)

#predico le probabilità per la sequenza ignota
prediction <- predict(neural_net, newdata=as.data.frame(raw_data[1:nrow(raw_data), 2:ncol(raw_data)], col.names=colnames(raw_data)[2:ncol(raw_data)], row.names=raw_data[1:nrow(raw_data), 1], stringsAsFactors=FALSE))

#stampo le probabilità predette
write.table(format(round(prediction, 4), , big.mark=",", scientific=FALSE), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
