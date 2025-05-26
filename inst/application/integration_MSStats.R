library(tidyverse)
library(prolfqua)
library(prophosqua)
library(dplyr)
library(stringr)
library(seqinr)

# integration of phospho-DIA
message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")


# parameters and thresholds
# variables
wu_id <- "CustomPhosphoAnalysis"
fgcz_project <- "p37875"
oid_fgcz <- "p37875"
fracti <- "Integration"
datetoday <- format(Sys.Date(), "%Y%m%d")
descri <- wu_id # "OneCondition"
(res_dir <- paste0("V3_", fgcz_project, "_", datetoday, "_", fracti, "_", descri))

path <- "."


# find path manually
# read back in results
tot_xlsx <- file.path(
  "DEA_proteome",
  "DEA_20250326_O37875_WU37875_complete_proteome_vsn",
  "Results_WU_37875_complete_proteome",
  "DE_WU37875_complete_proteome.xlsx"
)
phospho_xlsx <- "DEA_PTM/DEA_20250326_O37875_WU37875_enriched_vsn//Results_WU_37875_enriched/DE_WU37875_enriched.xlsx"


tot_res <- readxl::read_xlsx(path = tot_xlsx, sheet = "diff_exp_analysis")
# check how many proteins have no length ->  only rev proteins should have not length
tot_res$protein_length |>
  is.na() |>
  mean()
# these are the protein IDs without length! # should be only rev_ or proteins that are not found in this fasta
xx <- tot_res |>
  dplyr::filter(is.na(protein_length)) |>
  dplyr::pull(protein_Id)

rev_pattern <- "^rev_"
tot_res <- tot_res |> dplyr::filter(!grepl("FGCZCont", protein_Id))
tot_res <- tot_res |> dplyr::filter(!grepl("contam_", protein_Id))
tot_res <- tot_res |> dplyr::filter(!grepl(rev_pattern, protein_Id)) # here we take out all revs!
tot_res <- tot_res |> dplyr::filter(!contrast %in% c("ask13FCvsWTFC", "H6FCvsH1FC"))
ggplot(tot_res, aes(x = avgAbd)) +
  geom_histogram() +
  facet_wrap(~contrast)

phospho_res <- readxl::read_xlsx(path = phospho_xlsx, sheet = "diff_exp_analysis")
# drop contaminants and rev sequences
phospho_res <- phospho_res |> dplyr::filter(!grepl("FGCZCont", protein_Id))
phospho_res <- phospho_res |> dplyr::filter(!grepl("cont_", protein_Id))
phospho_res <- phospho_res |> dplyr::filter(!grepl(rev_pattern, protein_Id))
phospho_res <- phospho_res |> dplyr::filter(!contrast %in% c("ask13FCvsWTFC", "H6FCvsH1FC"))

# add here more phospho related things

# add things to phospho site
phospho_res2 <- phospho_res |> separate(site, c("proteinID", "SiteNlocation", "sequence"), remove = FALSE, sep = "_|~")
# needed downstream
phospho_res2 <- phospho_res2 |> mutate(posInProtein = as.numeric(gsub("[S|T|Y](\\d+)", "\\1", SiteNlocation)))
phospho_res2 <- phospho_res2 |> mutate(modAA = gsub(pattern = "([S|T|Y])\\d+", replacement = "\\1", SiteNlocation))
phospho_res2$modAA |> table()

phospho_res2$startModSite <- phospho_res2$posInProtein - 1
phospho_res2$endModSite <- phospho_res2$posInProtein + 1
# get the first modified AA

# num phos, numLoc and AllLocalized
phospho_res2$AllLocalized <- TRUE


# join total and phospho analysis
join_column <- c("fasta.id" = "protein_Id", "contrast", "description", "protein_length")
join_column %in% colnames(phosphoRes2)
join_column %in% colnames(totRes)

# do diff-diff and additional propagation of errors as suggested by MSstatsPTM
suffix_a <- ".site"
suffix_b <- ".protein"
combined_site_prot <- dplyr::left_join(phospho_res2, tot_res, by = join_column, suffix = c(suffix_a, suffix_b))

resDir
dir.create(resDir)

#
# write to excel
excel_result_list <- list()
excel_result_list$combinedStats <- combined_site_prot
nrow(combined_site_prot)
writexl::write_xlsx(excel_result_list, path = file.path(resDir, "Result_phosphoAndTotal.xlsx"))


# N-to-C plotting
# Function to determine the significance for plotting NtoC
# only proteins where in at least one contrast the protein is significant (fdrThreshold) are plotted
fdr_threshold <- 0.01
# build up candidate matrix for plotting
candidate_mat <- combined_site_prot[!is.na(combined_site_prot$FDR.site), ]
cand <- candidate_mat[candidate_mat$FDR.site < fdr_threshold, ]
# proteins with sites regulated in any of the contrasts.
sig_protein_hits <- unique(cand$protein_Id)

# look at one of the contrasts
combo_mat <- candidate_mat
combo_mat <- combo_mat |> dplyr::filter(protein_Id %in% sig_protein_hits)

# these columns are used later for N-to-C plots
combo_mat_min <- dplyr::select(
  combo_mat,
  c("protein_Id",
    "contrast",
    "protein_length",
    "site",
    "diff.protein",
    "diff.site",
    "FDR.site",
    "posInProtein",
    "startModSite",
    "endModSite",
    "AllLocalized",
    "modAA",
    model_site = "modelName.site"
  )
)

# nest all sites behind
combo_mat_min$site <- gsub("~.*", "", combo_mat_min$site)
combo_mat_min <- combo_mat_min |>
  dplyr::group_by(protein_Id, contrast, protein_length) |>
  tidyr::nest()
combo_mat_min$plot <- vector(mode = "list", nrow(combo_mat_min))

# fill all slots with plots

for (i in seq_len(nrow(combo_mat_min))) {
  print(i)
  combo_mat_min$plot[[i]] <- N_to_C_plot(
    combo_mat_min$data[[i]],
    combo_mat_min$protein_Id[[i]],
    combo_mat_min$protein_length[[i]],
    combo_mat_min$contrast[[i]]
  )
}



combo_mat_min <- ungroup(combo_mat_min)
combo_mat_min_wr <- select(combo_mat_min, contrast, plot, data)
tmp <- combo_mat_min_wr |>
  group_by(contrast) |>
  group_split()

# Do the plotting only for each contrast individually
if (TRUE) {
  for (j in seq_along(tmp)) {
    print(j)
    combo_mat_min <- tmp[[j]]
    pdf_fn <- paste0("Significant_Sites_Proteins_", combo_mat_min$contrast[j], "_NtoCplots.pdf")
    pdf(file.path(resDir, pdf_fn))
    for (i in seq_len(nrow(combo_mat_min))) {
      print(combo_mat_min$plot[[i]])
      grid::grid.newpage()
      table <- combo_mat_min$data[[i]] |> select(-all_of(c("startModSite", "endModSite", "AllLocalized")))
      table_grob <- gridExtra::tableGrob(table, theme = gridExtra::ttheme_default(base_size = 6))
      grid::grid.draw(table_grob)
    }
    dev.off()
  }
}


#### Integration data

# do diff-diff and additional propagation of errors as suggested by MSstatsPTM
combined_test_diff <- prophosqua::test_diff(phosphoRes2, totRes, join_column = join_column)

drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = fgcz_project,
  ORDERID = oid_fgcz,
  WORKUNITID = descri
)

prophosqua::copy_phospho_Integration()

# render html
rmarkdown::render("_Overview_PhosphoAndIntegration_site.Rmd",
  params = list(
    data = combined_test_diff,
    grp = drumm,
    phosres = phosphoRes
  ),
  output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE)
)
file.copy(
  from = "_Overview_PhosphoAndIntegration_site.html",
  to = file.path(resDir, "Result_phosphoAndIntegration.html")
)

# write to excel
excel_result_list <- list()
excel_result_list$combinedStats <- combined_test_diff
nrow(combined_test_diff)
writexl::write_xlsx(excel_result_list, path = file.path(resDir, "Result_phosphoAndTotalIntegration.xlsx"))


# Function to determine the significance for plotting NtoC
# only proteins where in at least one contrast the protein is significant (fdrThreshold) are plotted
fdr_threshold <- 0.01
# build up candidate matrix for plotting
candidate_mat <- combined_test_diff[!is.na(combined_test_diff$FDR_I), ]
cand <- candidate_mat[candidate_mat$FDR_I < fdr_threshold, ]

# proteins with sites regulated in any of the contrasts.
sig_protein_hits <- unique(cand$protein_Id)

# look at one of the contrasts
combo_mat <- candidate_mat
combo_mat <- combo_mat |> dplyr::filter(protein_Id %in% sig_protein_hits)
# how many true are ok -> position parsed properly -> this is irrelevant for site centric approach here
table(!is.na(combo_mat$posInProtein))
colnames(combo_mat)

# these columns are used later for N-to-C plots
combo_mat_min <- dplyr::select(
  combo_mat,
  c(
    "protein_Id",
    "contrast",
    "protein_length",
    "site",
    "diff_diff",
    "FDR_I",
    "posInProtein",
    "modAA"
  )
)
combo_mat_min$site <- gsub("~.*", "", combo_mat_min$site)
# nest all sites behind
combo_mat_min <- combo_mat_min |>
  dplyr::group_by(protein_Id, contrast, protein_length) |>
  tidyr::nest()
combo_mat_min$plot <- vector(mode = "list", nrow(combo_mat_min))
# check how many pages are plotted
length(unique(combo_mat_min$protein_Id))
combo_mat_min



# fill all slots with plots
for (i in seq_len(nrow(combo_mat_min))) {
  combo_mat_min$plot[[i]] <- N_to_C_plot_integrated(
    combo_mat_min$data[[i]],
    combo_mat_min$protein_Id[[i]],
    combo_mat_min$protein_length[[i]],
    combo_mat_min$contrast[[i]]
  )
}

combo_mat_min_wr <- select(ungroup(combo_mat_min), contrast, plot, data)
tmp <- combo_mat_min_wr |>
  group_by(contrast) |>
  group_split()



tmp <- combo_mat_min_wr |>
  group_by(contrast) |>
  group_split()
combo_mat_min <- tmp[[1]]
combo_mat_min$data[[2]]
# Do the plotting only for each contrast individually
if (TRUE) {
  for (j in seq_along(tmp)) {
    print(j)
    combo_mat_min <- tmp[[j]]
    pdf_fn <- paste0("SignificantSites_Normalized", unique(combo_mat_min$contrast), "_NtoCplots.pdf")
    pdf(file.path(resDir, pdf_fn))
    for (i in seq_len(nrow(combo_mat_min))) {
      print(i)
      print(combo_mat_min$plot[[i]])
      grid::grid.newpage()
      table <- combo_mat_min$data[[i]]
      table_grob <- gridExtra::tableGrob(table, theme = gridExtra::ttheme_default(base_size = 6))
      grid::grid.draw(table_grob)
    }
    dev.off()
  }
}
