# Utility function to read multiple identical tables and concatenate them.
read.identical <- function (file.names, header.columns, data.columns, file.labels = basename(file.names), col.names=NULL) 
{
    results = NULL
    for (i in 1:length(file.names)) {
        file.name = file.names[i]
        file.label = file.labels[i]
        file.data = read.table(file.name, sep = "\t", header = FALSE, 
            stringsAsFactors = FALSE)
        
        if(!is.null(col.names)) {
            colnames(file.data) = col.names
        }
        
        if (is.null(results)) {
            results = file.data[, header.columns, drop=FALSE]
        }
        colnames(file.data) <- paste(file.label, colnames(file.data), 
            sep = ".")
        results = cbind(results, file.data[, data.columns, drop=FALSE])
    }
    return(results)
}

# Figure out input and output paths.
base.path = commandArgs()[1]
count.path = file.path(base.path, "counts")
input.files = list.files(count.path)

# Read all counts, get them into a matrix without header columns.
all.counts = read.identical(file.path(base.path, "counts", input.files), 1, 2, gsub(".counts", "", input.files), col.names=c("ENSEMBL", "Count"))
rownames(all.counts) = all.counts[,1]
all.counts <- all.counts[,-1]

# Count how many reads were assigned to a transcript.
assigned = apply(all.counts[grepl("ENSG", rownames(all.counts)),], 2, sum)
all.counts["__assigned",] = assigned

# Gather all non-transcript rows (__*)
summary.counts = all.counts[!grepl("ENSG", rownames(all.counts)),]

# Generate a summary table.
total.reads = apply(summary.counts, 2, sum)
summary.table = data.frame(Total = total.reads,
                           Mapped = as.numeric(total.reads - summary.counts["__not_aligned",]),
                           NoFeature = as.numeric(summary.counts["__no_feature",]),
                           MultipleAlignments = as.numeric(summary.counts["__alignment_not_unique",]),
                           Assigned = as.numeric(summary.counts["__assigned",]))
write.table(summary.table, file=file.path(base.path, "SummaryTable.txt"), sep="\t", row.names=FALSE, col.names=TRUE)