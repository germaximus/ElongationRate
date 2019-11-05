library(data.table)
library(magrittr)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#------------------------------------------ Define some useful functions -----------------------------------------------------------------------
linkage <- function(gff) { # creates a 2-column table with children->parent linkages. Takes original gff annotation as its argument.
  output <- apply(gff[,9], 1, function(x) {
    id <- substr(x, 4, regexpr(';', x) - 1)
    if(regexpr('Parent=', x)[[1]] > 0) { parent <- substr(x, regexpr('Parent=', x) + 7, gregexpr(';', x)[[1]][2] -1); return(c(id, parent))}
    else {return(c(id, "Primary"))}
  }) %>% t() %>% as.data.frame(., stringsAsFactors = F) %>% setNames(., c("ID", "Parent1"))
  return(output)
}

chain <- function(x) { # reconstitutes a full chain of parents for every child record. Takes the output of the linkage function as its argument.
  temp <- x[match(x[, ncol(x)], x$ID), "Parent1"]
  temp[is.na(temp)] <- "Primary"
  temp <- as.data.frame(temp, stringsAsFactors = F) %>% setNames(., paste0("Parent", ncol(x)))
  print(length(unique(temp[,1])))
  if(length(unique(temp[,1])) == 1) {return(x)}
  else{ return(chain(cbind(x, temp))) }
}

remove.features <- function(gff, type, parents){ # removes unwanted feature types from the gff file together with their children and parents
              id <- gff[unlist(gff[,3]) %in% type, 9] %>% apply(., 1, function(x) {substr(x, 4, regexpr(';', x) - 1)})
              #parents <- chain(linkage(gff))
              primary_id <- apply(parents[match(id, parents$ID), ], 1, function(x) {
                record <- x[x != "Primary"] 
                record <- record[length(record)]
                return(record)
              }) %>% unname() %>% unique()
              
              discard_lines <- sapply(parents, function(x) {  x %in% primary_id    }) %>% apply(., 1, function(x) { any(x) }) 
              output <- gff[!discard_lines, ]
              return(output)
}

extract.id <- function(gff, type, level = "Primary", parents){ # simplified version of remove.features() that only reports top-level (Primary) or their own local IDs for given feature types
              id <- gff[unlist(gff[,3]) %in% type, 9] %>% apply(., 1, function(x) {substr(x, 4, regexpr(';', x) - 1)})
              if(level != "Primary") {return(id)}
              else { #parents <- chain(linkage(gff))
                     primary_id <- apply(parents[match(id, parents$ID), ], 1, function(x) {
                        record <- x[x != "Primary"] 
                        record <- record[length(record)]
                        return(record)
                      }) %>% unname() %>% unique()
              return(primary_id)
              }
}

#----------------------------------------------------------------------------------------------------------------------------------------
# Load GFF annotation file
gff    <- fread(file="GRCm38.p6.Refseq.gff", skip = 9, stringsAsFactors = F, header = F, fill = T, na.strings = c("", "NA"), sep="\t") %>% na.omit() #deals with unwanted #comment lines in the gff
gff$ID <- apply(gff[,9], 1, function(x) { id <- substr(x, 4, regexpr(';', x) - 1) })

# create a table of parents-children relations
parents.table <- chain(linkage(gff))

# Create a list of top parents with all childrens listed (Caution: takes ~1 hour) # Save the R object for future use to avoid re-calculations.
parents.tree <- lapply(unique(parents.table[parents.table$Parent1 == 'Primary','ID']), function(x) { return(NULL) }) %>% setNames(unique(parents.table[parents.table$Parent1 == 'Primary', 'ID']))
for(i in 1:nrow(parents.table)) {
          if (!isTRUE(parents.table[i,'ID'] %in% names(parents.tree)) ) {
            children <- parents.table[i, -c(1)] %>% unname() %>% unlist()
            parent   <- children[children %in% names(parents.tree)]
            parents.tree[[parent]] <- c(parents.tree[[parent]], parents.table[i, 'ID']) %>% unique()   }
} 

#saveRDS(parents.tree, file = "parents_tree.rds")
#parents.tree <- readRDS('parents_tree.rds')

# fix gene boundaries to be the same as the span of children features. The discrepancy occured when I removed Gnomon records. Some gene names were shared between BestRefseq and Gnomon as 'BestRefSeq%2CGnomon'. Their boundaries are wider than corresponding BestRefSeq childs.
setindex(gff, ID)
for(name in names(parents.tree)) {
    slice <- gff[name, on = 'ID']

    if(nrow(slice) == 1){
      gene_start <- slice[, V4] 
      gene_end   <- slice[, V5] 

      if(!is.null(parents.tree[[name]]))    {
        children_start <- gff[parents.tree[[name]], V4, on = 'ID'] %>% min()
        children_end   <- gff[parents.tree[[name]], V5, on = 'ID'] %>% max()
        
        if(children_start > gene_start) { gene_start = children_start }
        if(children_end   < gene_end  ) { gene_end   = children_end   }
      } 
      gff[name, on = 'ID', V4 := gene_start]
      gff[name, on = 'ID', V5 := gene_end]
    }
} 
         #con <- file("GRCm38.p6.Refseq.gff", "r")
         #header <- readLines(con, n = 7)
         #write.table(header, file = "GRCm38.p6.Refseq.fixed.gff", col.names = F, row.names = F, quote = F)
         #write.table(gff[,1:9], file = "GRCm38.p6.Refseq.fixed.gff", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
         #close(con); rm(con)

# Remove non-coding features
# check all present top level features: table(gff[,3])
# Tips: NCBI RefSeq. Difference between 'mRNA' and 'transcript' features is the former have CDS and the latter do not.
# Gene --> mRNA --> CDS and exons.
# Gene --> transcript --> exons
# Therefore, do not remove 'transcripts' otherwise you remove their gene parents.

gff2 <- remove.features(gff, c('antisense_RNA','lnc_RNA','region','snRNA','cDNA_match','match','RNase_MRP_RNA','telomerase_RNA','miRNA','RNase_P_RNA','primary_transcript','rRNA','tRNA','scRNA','Y_RNA','guide_RNA','pseudogene','snoRNA'), parents.table)
con <- file("GRCm38.p6.Refseq.gff", "r")
header <- readLines(con, n = 7)
write.table(header, file = "GRCm38.p6.Refseq.coding.gff", col.names = F, row.names = F, quote = F)
write.table(gff2[,1:9], file = "GRCm38.p6.Refseq.coding.gff", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
close(con); rm(con)































