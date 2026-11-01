library(msigdbr)
r <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "C2",
                                          subcategory = "CP:REACTOME"))
want <- c(IFNA = "REGULATION_OF_IFNA", IL4_IL13 = "INTERLEUKIN_4_AND_INTERLEUKIN_13",
          PLATELET_COLLAGEN = "PLATELET_ADHESION_TO_EXPOSED_COLLAGEN",
          SARS = "POTENTIAL_THERAPEUTICS_FOR_SARS", TRYPTOPHAN = "TRYPTOPHAN_CATABOLISM",
          UPR = "UNFOLDED_PROTEIN_RESPONSE", ASP_ASN = "ASPARTATE_AND_ASPARAGINE",
          BCAA = "BRANCHED_CHAIN_AMINO_ACID_CATABOLISM",
          FZD = "REGULATION_OF_FZD_BY_UBIQUITINATION",
          SULFUR = "SULFUR_AMINO_ACID_METABOLISM", PEROX = "PEROXISOMAL_PROTEIN_IMPORT")
nm <- unique(r$gs_name)
for (k in names(want)) {
  hit <- grep(want[[k]], nm, value = TRUE)
  cat(sprintf("%-18s -> %s\n", k, paste(hit, collapse = " | ")))
  if (length(hit) == 1)
    writeLines(unique(r$gene_symbol[r$gs_name == hit]),
               file.path("data/signatures", paste0("Linggi_", k, ".txt")))
}
list.files("data/signatures")
fix <- c(PLATELET_COLLAGEN = "REACTOME_PLATELET_ADHESION_TO_EXPOSED_COLLAGEN",
         UPR               = "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
         BCAA              = "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM")
for (k in names(fix))
  writeLines(unique(r$gene_symbol[r$gs_name == fix[[k]]]),
             file.path("data/signatures", paste0("Linggi_", k, ".txt")))
length(list.files("data/signatures"))   # should print 11
source("06_signatures_interaction.R")
