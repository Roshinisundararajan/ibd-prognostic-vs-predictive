install.packages("msigdbr")
library(msigdbr)
h <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
sets <- split(h$gene_symbol, h$gs_name)
dir.create("data/genesets", showWarnings = FALSE)
writeLines(vapply(names(sets), function(n)
  paste(c(n, "MSigDB_Hallmark", unique(sets[[n]])), collapse = "\t"),
  character(1)),
  "data/genesets/hallmark.gmt")
length(sets)                 # should print 50
packageVersion("msigdbr")    # note this for your methods
