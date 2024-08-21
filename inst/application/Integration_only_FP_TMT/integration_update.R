message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")

# annotation and descrison based on:
# https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


resDir <- "p35593_uniprot_paired/"
totXlsx <- "p35593_uniprot_paired/DEA_20240821_OI_35593_WU_WholeProtUniprot_paired_vsn/Results_DEA_WUWholeProtUniprot_paired/DE_Groups_vs_Controls_WUWholeProtUniprot_paired.xlsx"
phosXlsx <- "p35593_uniprot_paired/DEA_20240821_OI_35593_WU_EnrichedPhosUniprot_paired_vsn/Results_DEA_WUEnrichedPhosUniprot_paired/DE_Groups_vs_Controls_WUEnrichedPhosUniprot_paired.xlsx"

totRes <- readxl::read_xlsx(path = totXlsx, sheet = "diff_exp_analysis")
totRes$protein_length |> is.na() |> mean()
totRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id)

totRes <- totRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl("contam_", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl("rev_", protein_Id))

rev_pattern = "^rev_"

phosRes <- readxl::read_xlsx(path = phosXlsx, sheet = "diff_exp_analysis")
phosRes$protein_length |> is.na() |> mean()


# drop contaminants and rev sequences
phosRes <- phosRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
phosRes <- phosRes |> dplyr::filter(!grepl("cont_", protein_Id))
phosRes <- phosRes |> dplyr::filter(!grepl("rev_", protein_Id))

phosRes <- phosRes |> dplyr::mutate(
  AllLocalized = (NumPhos == LocalizedNumPhos)
)

phosRes$AllLocalized |> mean()

fasta_file = "o35593_prot_ionquant/2024-07-03-decoys-contam-UP000186698.fasta"

get_sequence_windows <- function(phosRes, fasta_file, rev_pattern = "rev_") {
  uniqueProtPepSeq <- phosRes |> dplyr::filter(AllLocalized == TRUE) |>
    dplyr::select(protein_Id, site, PhosSites) |> dplyr::distinct()

  uniqueProtPepSeq <- uniqueProtPepSeq |> tidyr::separate_longer_delim(PhosSites, delim = ";")
  uniqueProtPepSeq <- uniqueProtPepSeq |> dplyr::mutate(posInProtein = as.integer(stringr::str_remove(PhosSites, "^[A-Z]")))
  uniqueProtPepSeq <- uniqueProtPepSeq |> dplyr::mutate(AA = stringr::str_remove(PhosSites, "\\d+"))


  fasta <- prolfquapp::get_annot_from_fasta(fasta_file, rev = rev_pattern, include_seq = TRUE)

  uniqueProtPepSeq <- dplyr::inner_join(uniqueProtPepSeq, fasta, by = c(protein_Id = "proteinname") )
  window_size <- 15
  uniqueProtPepSeq <- uniqueProtPepSeq |>
    dplyr::mutate(
      padded_sequence = paste0(strrep("X", window_size), sequence, strrep("X", window_size)),
      posInPaddedSeq = posInProtein + window_size,  # Adjust position due to padding
      posStart = posInPaddedSeq - window_size,
      posEnd = posInPaddedSeq + window_size,
      sequence_window = substr(padded_sequence, start = posStart, stop = posEnd)
    ) |>
    dplyr::select(-padded_sequence, -posInPaddedSeq, -posStart, -posEnd)

  uniqueProtPepSeq$sequence <- NULL
  uniqueProtPepSeq
}

seq_window <- get_sequence_windows(phosRes, fasta_file, rev_pattern)

join_column <- c("fasta.id" = "protein_Id", "contrast", "description", "protein_length", "nr_tryptic_peptides")


reverse_join_column <- function(join_column){
  reverse_join_column <- vector(mode = "character", length(join_column))
  for (i in seq_along(join_column)) {
    reverse_join_column[i] <- if (names(join_column)[i] != "") {  names(join_column)[i]} else { join_column[i]}
    names(reverse_join_column)[i] <- if (names(join_column)[i] != "") { join_column[i]} else {""}
  }
  return(reverse_join_column)
}


test_diff <- function(phosRes, totRes, join_column = c("protein_Id", "contrast","description", "protein_length", "nr_tryptic_peptides")){
  test_diff <- prophosqua::test_diff_diff(phosRes,totRes, by = join_column)
  test_diff$measured_In <- "both"

  removed_from_site <- dplyr::anti_join(phosRes, totRes, by = join_column )
  removed_from_site$measured_In <- "site"

  removed_from_prot <- dplyr::anti_join(totRes, phosRes, by = reverse_join_column(join_column))
  removed_from_prot$measured_In <- "prot"

  common_columns <- setdiff(intersect(colnames(removed_from_site), colnames(removed_from_prot)),c(join_column, "measured_In"))
  removed_from_site_renamed <- removed_from_site |>
    dplyr::rename_with(~ paste0(., ".site"), tidyselect::all_of(common_columns))
  removed_from_prot_renamed <- removed_from_prot |>
    dplyr::rename_with(~ paste0(., ".protein"), dplyr::all_of(common_columns))

  combined_test_diff <- dplyr::bind_rows(test_diff , removed_from_site_renamed , removed_from_prot_renamed)
  return(combined_test_diff)
}

combined_test_diff <- test_diff(phosRes, totRes, join_column = c("fasta.id" = "protein_Id", "contrast","description", "protein_length", "nr_tryptic_peptides"))
combined_test_diff$AllLocalized |> mean(na.rm = TRUE)

phosRes <- phosRes |> dplyr::mutate(AA = substr(PhosSites, 1, 1))


drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = "p35593",
  ORDERID = "p35593",
  WORKUNITID = "Xenbase Database")
rmarkdown::render("_Overview_PhosphoAndIntegration_WEW.Rmd", params = list(data = combined_test_diff, grp = drumm, phosres = phosRes), output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE))
file.copy(from = "_Overview_PhosphoAndIntegration_WEW.html", to = file.path(resDir, "PhosphoAndIntegration.html"))

# separate for each constrast

# write to excel
excelResultList <- list()
excelResultList$combinedStats <- combined_test_diff
excelResultList$seq_window <- seq_window
writexl::write_xlsx(excelResultList, path = file.path(resDir,"PhosphoAndIntegration.xlsx"))
#

# Function to determine the significance annotation
fdrThreshold = 0.2


candidateMat <- combined_test_diff[!is.na(combined_test_diff$FDR.site), ]
cand <- candidateMat[candidateMat$FDR.site < fdrThreshold ,]
# proteins with sites regulated in any of the contrasts.
mySigProteinHits <- unique(cand$protein_Id)

# look at one of the contrasts
comboMat <- candidateMat
comboMat <- comboMat |> dplyr::filter(protein_Id %in% mySigProteinHits)
candidateMat$AllLocalized |> mean()

# add data frame with sequence windows
comboMat <- dplyr::inner_join(comboMat, seq_window)


# select columns necessary for plotting
comboMat_min <- dplyr::select(
  comboMat,
  c("protein_Id",
    "contrast",
    "protein_length",
    "site",
    "diff.protein",
    "diff.site",
    "FDR.site",
    "modAA" = "AA",
    "posInProtein",
    "startModSite",
    "endModSite",
    "AllLocalized",
    model_site = "modelName.site"
  ))


comboMat_min <- comboMat_min |> dplyr::group_by(protein_Id, contrast, protein_length) |> tidyr::nest()
comboMat_min$plot <- vector(mode = "list", nrow(comboMat_min))

library(tidyverse)

for (i in 1:nrow(comboMat_min)) {
  print(i)
  comboMat_min$plot[[i]] <- prophosqua::N_to_C_plot(comboMat_min$data[[i]],
                                                    comboMat_min$protein_Id[[i]],
                                                    comboMat_min$protein_length[[i]],
                                                    comboMat_min$contrast[[i]])
}

head(comboMat_min)


pdf(file.path(resDir, "NtoCplots2.pdf"))

for (i in 1:nrow(comboMat_min)) {
  print(comboMat_min$plot[[i]])
  grid::grid.newpage()
  table <- comboMat_min$data[[i]]
  table <- table |> select(-all_of(c( "startModSite", "endModSite", "AllLocalized")))
  table_grob <- gridExtra::tableGrob(table, theme = gridExtra::ttheme_default(base_size=6))
  grid::grid.draw(table_grob)
}
dev.off()
