#' read AA fasta file, identify decoy entries and return fw fasta object only
#' @param FastaFileName filename of fasta file
#' @param decoyPattern pattern for decoy accessions
#' @return fastaObject without decoy sequences
#' @export
#'
readDecoyFastaNreturnFwOnly <- function(FastaFileName, decoyPattern = "REV_") {
  myFasta_decoy <- seqinr::read.fasta(file = FastaFileName, seqtype = "AA", as.string = TRUE)
  seqNames <- getName(myFasta_decoy)
  idx_rev <- which(str_count(string = seqNames, pattern = decoyPattern)>0)
  nodeoySeq <- myFasta_decoy[-idx_rev]
  return(nodeoySeq)
}



#' prepare longformat from FragPipe Site Centric Phospho output called STY_79.9663.tsv
#' @param poi protein of interest from FP-site centric output
#' @param soi site of interest (parsed from index in FP-site centric output)
#' @param fastaObject read by seqinR
#' @param seqWindowOffset Sequence window offset plus minus from site
#' @return string of X characters with p-Site in the middle
#' @export
#'
getSequenceWindowForLogo <- function(poi, soi, fastaObject = myFasta, seqWindowOffset = 15) {
  # find poi in fasta -> here we do not handle cases where it might appear multiple times in the full accession
  f_id <- which(str_count(string = getName(fastaObject), pattern = poi)>0)
  if((soi-seqWindowOffset) < 0) {
    soiStart <- 0
    } else {
    soiStart <- soi-seqWindowOffset
  }
  if((soi+seqWindowOffset) > nchar(fastaObject[[f_id]][1])) {
    soiEnd <- nchar(fastaObject[[f_id]][1])
    } else {
    soiEnd <- soi+seqWindowOffset
  }
  mySeq <- as.character(getSequence(getFrag(object = fastaObject[[f_id]], begin = soiStart, end = soiEnd, as.string=TRUE), as.string = TRUE))
  return(mySeq)
}



#' prepare longformat from FragPipe Site Centric Phospho output called STY_79.9663.tsv
#' @param styphosphopeptideTable original STY_79.9663 file from FragPipe
#' @param locProbThreshold Threshold on best Localization Probability
#' @param locPnames columns for individual  localization site probability
#' @param protIDcol id column (index = proteinName_site)
#' @param intnames tag to identify intensity columns
#' @param maxlfqnames tag to identify MaxLFQ intensity columns
#' @return data.frame in lonng  format
#' @export
#'
tidy_FragPipe_phospho_STYfile <- function (styphosphopeptideTable, as_list = FALSE,
                                           locPnames = c(" Localization Probability"),
                                           locProbThreshold = 0.75,
                                           protIDcol = "Index", intnames = c(" Intensity"),
                                           maxlfqnames = c("MaxLFQ Intensity"))
{
  if (is.character(styphosphopeptideTable) && file.exists(styphosphopeptideTable)) {
    Cprotein <- as_tibble(read.csv(styphosphopeptideTable, header = TRUE,
                                   sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
  }
  else if ("tbl_df" %in% class(styphosphopeptideTable)) {
    Cprotein <- styphosphopeptideTable
  }
  else {
    stop(class(styphosphopeptideTable), " not supported.")
  }
  #cnam <- gsub("Total Razor ", "Total ", gsub("Unique Razor ",
  #                                            "Unique ", gsub(" Intensity$", " Razor Intensity", gsub(" Spectral Count$",
  #                                                                                                    " Razor Spectral Count", colnames(Cprotein)))))


  # replace "MaxLFQ Intensity" ->"MaxLFQIntensity" for parsing
  maxlfqnames <- gsub(x = maxlfqnames, pattern = "MaxLFQ Intensity", replacement = "MaxLFQIntensity")
  cnam <- gsub(x = colnames(Cprotein), pattern = "MaxLFQ Intensity", replacement = "MaxLFQIntensity")
  colnames(Cprotein) <- cnam

  # Filter here for Best Localization Probability and then leave it away completely
  Cprotein <- Cprotein |> dplyr::filter(Cprotein$`Best Localization Probability` > locProbThreshold)
  cnam <- cnam[1:which(cnam == "Best Localization Probability")]
  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- dplyr::select(Cprotein, all_of(cnam))

  # jg: not clear why witold did it this way, issue with Best Loc Prob!
  extractDataLong <- function(Cprotein, what = "Intensity",
                              butNot = NULL) {
    cols <- colnames(Cprotein)
    cols <- setdiff(grep(paste0(what, "$"), cols, value = TRUE),
                    if (is.null(butNot)) {
                      NULL
                    }
                    else {
                      grep(butNot, cols, value = TRUE)
                    })
    gg <- dplyr::select(Cprotein, all_of(protIDcol), all_of(cols))
    gg <- tidyr::pivot_longer(gg, cols = dplyr::ends_with(what),
                              names_to = "raw.file", values_to = what)
    # gg <- dplyr::mutate(gg, raw.file = gsub(paste0(".", what,
    gg <- dplyr::mutate(gg, raw.file = gsub(paste0("", what,
                                                   "$"), "", .data$raw.file))
    gg
  }
  res <- vector(mode = "list", length = length(c(intnames)))
  names(res) <- c(intnames)
  for (i in seq_along(c(intnames))) {
    message("DD: ", c(intnames)[i])
    res[[c(intnames)[i]]] <- extractDataLong(Cprotein,
                                                       what = c(intnames)[i], butNot = "maxlfq")
  }
  if (sum(grepl(".MaxLFQ.", colnames(Cprotein))) > 0) {
    res_maxlfq <- vector(mode = "list", length(maxlfqnames))
    names(res_maxlfq) <- maxlfqnames
    for (i in seq_along(maxlfqnames)) {
      message("DD: ", maxlfqnames[i])
      res_maxlfq[[maxlfqnames[i]]] <- extractDataLong(Cprotein,
                                                      what = maxlfqnames[i], butNot = NULL)
    }
    res <- c(res, res_maxlfq)
  }
  if (as_list) {
    return(res)
  }
  # sql_inner_join <- function(x, y) {
  #   inner_join(x, y, multiple = "all")
  # }
  #merged <- Reduce(sql_inner_join, res)
  #merged <- inner_join(annot, merged, multiple = "all")
  merged <- left_join(x = res$` Intensity`, y = annot)
  # fix raw-file names for maxLFQInt because we need the space for proper Intensity parsing
  res$MaxLFQIntensity$raw.file <- gsub(x = res$MaxLFQIntensity$raw.file, pattern = " ", replacement = "")
  merged <- left_join(x = merged, y = res$MaxLFQIntensity)
  colnames(merged) <- gsub(x = colnames(merged), pattern = " ", replacement = "")
  colnames(merged) <- tolower(make.names(colnames(merged)))
  return(merged)
}





# Do the NtoCplot
# 2023-07-05 adapted for FP outputs
# this undocumented function is used to plot the phospho peptides from NtoC-term with indicated log2-fold-change-bars
# the function is used inside another documented function -> the 2019 function should be used from 2019 onwards
.plotProteinNPhosphoPeptidesNtoC_FragPipe <- function(protName, globalProtFC, POI_tableProtein, POI_tablePeptide,
                                                  protColor = "blue", protWtdh = 4, maxX = 1200, protNCoffset = 40,
                                                  pepSigStarOffset = 0.05, pepSigStarSize = 2, fdrRelaxedThreshold = 0.2,
                                                  fdrStringentThreshold = 0.05, greyLineLog2Threshold = 1) {
  myYaxisLimiter <- max(ceiling(max(abs(POI_tablePeptide$diff.x))), ceiling(max(abs(POI_tableProtein$diff.y))), na.rm = TRUE)
  miny <- -myYaxisLimiter
  maxy <- myYaxisLimiter
  # get positions
  Phos_PositionsInProteins <- vector(length=nrow(POI_tablePeptide))
  for (i in 1:nrow(POI_tablePeptide)) {
    Phos_PositionsInProteins[i] <-  as.numeric(POI_tablePeptide$posInProtein[i])
   }
  scaleX <- 1/max(Phos_PositionsInProteins)*1200
  POI_tablePeptide$PositionsINproteins <- Phos_PositionsInProteins
  phosphoPlotTitle <- paste("Prot: ",protName,"\n Peptides from: ", POI_tablePeptide$proteinid[1],"# phospho = ",
                            length(POI_tablePeptide$diff.x), sep=" ")
  # Do the plot start with white
  plot(c(-100,0,1400, 0,0), c(globalProtFC,0,0,miny, maxy), pch=".", col="white", main=phosphoPlotTitle, ylab="log2FC",
       xlab="Full length protein", xaxt='n')
  axis(side = 1,at = c(-100,0,1200),labels = c("protFC","0","100%"))
  #grey box for prot FC
  rect(xleft = -150,ybottom = -100,xright = -50,ytop = 100, col = "lightgrey",border = FALSE)
  if (is.na(globalProtFC)) {
    points(-100,0,pch="0", col="orange", cex=3)
  } else {
    segments(-100,0,-100, globalProtFC, col=protColor, lwd=protWtdh)
  }

  # set the stage for the protein
  segments(0,0,maxX,0, col="black", lwd=1.5)
  points(-protNCoffset,0,pch="N")
  points(maxX+protNCoffset,0, pch="C")

  # plot all peptides
  #make modAA as factor
  POI_tablePeptide$modAA <- as.factor(POI_tablePeptide$AA)

  #make modelName a factor for lty -> strange error
  POI_tablePeptide$modelName.x <- as.factor(POI_tablePeptide$modelName.x)

  for (i in 1:length(POI_tablePeptide$PositionsINproteins)) {
    # plot lines
    segments(POI_tablePeptide$PositionsINproteins[i]*scaleX, 0,
             POI_tablePeptide$PositionsINproteins[i]*scaleX, POI_tablePeptide$diff.x[i],
             col=as.numeric(POI_tablePeptide$modAA[i]), lwd=1, lty=abs(as.numeric(POI_tablePeptide$modelName.x[i])-2)+1) # quite a hack for lty to get 2 = 1

    # decoration!
    # site is stringent significant
    if (POI_tablePeptide$FDR.x[i] < fdrStringentThreshold) {
            if(POI_tablePeptide$diff.x[i] > 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                        scaleX, POI_tablePeptide$diff.x[i]+pepSigStarOffset, pch="*",
                                                      cex=pepSigStarSize,
                                                      col=as.numeric(POI_tablePeptide$modAA[i]))

      if(POI_tablePeptide$diff.x[i] < 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                  scaleX, POI_tablePeptide$diff.x[i]-pepSigStarOffset, pch="*",
                                                cex=pepSigStarSize,
                                                col=as.numeric(POI_tablePeptide$modAA[i]))
    }
    else if (POI_tablePeptide$FDR.x[i] < fdrRelaxedThreshold) {
      if(POI_tablePeptide$diff.x[i] > 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                  scaleX, POI_tablePeptide$diff.x[i]+pepSigStarOffset, pch="+",
                                                cex=pepSigStarSize,
                                                col=as.numeric(POI_tablePeptide$modAA[i]))

      if(POI_tablePeptide$diff.x[i] < 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                  scaleX, POI_tablePeptide$diff.x[i]-pepSigStarOffset, pch="+",
                                                cex=pepSigStarSize,
                                                col=as.numeric(POI_tablePeptide$modAA[i]))
    }
    else if (POI_tablePeptide$FDR.x[i] > fdrRelaxedThreshold) {
      if(POI_tablePeptide$diff.x[i] > 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                  scaleX, POI_tablePeptide$diff.x[i]+pepSigStarOffset, pch="x",
                                                cex=pepSigStarSize/3,
                                                col=as.numeric(POI_tablePeptide$modAA[i]))

      if(POI_tablePeptide$diff.x[i] < 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                  scaleX, POI_tablePeptide$diff.x[i]-pepSigStarOffset, pch="x",
                                                cex=pepSigStarSize/3,
                                                col=as.numeric(POI_tablePeptide$modAA[i]))
      if(POI_tablePeptide$diff.x[i] == 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                  scaleX, POI_tablePeptide$diff.x[i], pch="+",
                                                cex=pepSigStarSize/3,
                                                col=as.numeric(POI_tablePeptide$modAA[i]))
      }

    }
  legend("bottomright", legend = unique(POI_tablePeptide$modAA),
         col=unique(as.numeric(POI_tablePeptide$modAA)),
         text.col = unique(as.numeric(POI_tablePeptide$modAA)), lwd=1)
  legend("topright", legend = c("dashed line -> imputed model", "below stringent FDR threshold", "below relaxed FDR threshold", "above FDR threshold"),
         col = c("orange", "black", "black", "black"), pch=c(" ","*", "+", "x"), cex = 0.7, text.col = c("orange", "black", "black", "black"))
  abline(h=greyLineLog2Threshold, col="grey", lty=2)
  abline(h=-greyLineLog2Threshold, col="grey", lty=2)
}



# 2023-07-05 adapted for FP outputs
#' Generate PDFs with NtoCplots for all candidates (using all phospho peptides even if not significant)
#' @param globalNphosphoCombinedNResultMatrix combinded result matrix with phospho and total protein result
#' @param candidateMatrix contains the significant phospho peptides for which all proteins should be drawn
#' @param expName this name is used to label the pdf file that is generated with this function
#' @export
#'
generateNtoCProteinPDFsWithPhosphoPeptides_FragPipe <- function(globalNphosphoCombinedNResultMatrix, candidateMatrix, expName) {
  # some static values to fill NAs
  # linValueFiller <- 1
  # qModForOneSiderCandidates <- 0.001
  pdfFileName <- paste("SignificantProtein_NtoCplot_", expName, ".pdf", sep="")
  mySigProteinHits <- rle(as.vector(sort(candidateMatrix$proteinid)))
  # Do open pdf here
  pdf(pdfFileName,10,10)
  for (j in 1:length(mySigProteinHits$values)) {
    rm(POI_matrix)
    POI <- mySigProteinHits$values[j]
    # Extract Relevant Lines from Full table again
    POI_matrix <- globalNphosphoCombinedNResultMatrix[which(globalNphosphoCombinedNResultMatrix$proteinid == POI), ]
    if(is.na(POI_matrix$p.value.y[1])) {
      message(paste("Protein Globally NOT quantified. Working on: ", unique(POI_matrix$proteinid)))
      MyProteinName_poi <- "Protein_globally_NOT_quantified"
      POI_protFC <- NA
    } else {
      message(paste("We got it globally. Working on: ", unique(POI_matrix$proteinid)))
      POI_protFC <- mean(POI_matrix$diff.y)
      MyProteinName_poi <- unique(POI_matrix$proteinid)
    }

        # Do the plot
    # split matrix in protein part and peptide part
    idx_split <- grep(x = colnames(POI_matrix), pattern = "description.y")
    POI_phosPeps<- POI_matrix[,1:(idx_split-1)]
    POI_globProts<- POI_matrix[,(idx_split-1):ncol(POI_matrix)]
    # Do the plotting
    .plotProteinNPhosphoPeptidesNtoC_FragPipe(protName = MyProteinName_poi, globalProtFC = POI_protFC,
                                          POI_tableProtein = POI_globProts, POI_tablePeptide = POI_phosPeps)
  }
  dev.off()
}


# 2023-07-05 adapted for FP outputs
#' Generate PDFs with NtoCplots for all candidates (using all phospho peptides even if not significant)
#' @param globalNphosphoCombinedNResultMatrix combinded result matrix with phospho and total protein result
#' @param candidateMatrix contains the significant phospho peptides for which all proteins should be drawn
#' @param expName this name is used to label the pdf file that is generated with this function
#' @export
#'
generateNtoCProteinPDFsWithPhosphoPeptides_FragPipeTMT <- function(globalNphosphoCombinedNResultMatrix, candidateMatrix, expName) {
  # some static values to fill NAs
  # linValueFiller <- 1
  # qModForOneSiderCandidates <- 0.001
  pdfFileName <- paste("SignificantProtein_NtoCplot_", expName, ".pdf", sep="")
  mySigProteinHits <- rle(as.vector(sort(candidateMatrix$protein_Id)))
  # Do open pdf here
  pdf(pdfFileName,10,10)
  for (j in 1:length(mySigProteinHits$values)) {
    rm(POI_matrix)
    POI <- mySigProteinHits$values[j]
    # Extract Relevant Lines from Full table again
    POI_matrix <- globalNphosphoCombinedNResultMatrix[which(globalNphosphoCombinedNResultMatrix$protein_Id == POI), ]
    if(is.na(POI_matrix$p.value.y[1])) {
      message(paste("Protein Globally NOT quantified. Working on: ", unique(POI_matrix$protein_Id)))
      MyProteinName_poi <- "Protein_globally_NOT_quantified"
      POI_protFC <- NA
    } else {
      message(paste("We got it globally. Working on: ", unique(POI_matrix$protein_Id)))
      POI_protFC <- mean(POI_matrix$diff.y)
      MyProteinName_poi <- unique(POI_matrix$protein_Id)
    }

    # Do the plot
    # split matrix in protein part and peptide part
    idx_split <- grep(x = colnames(POI_matrix), pattern = "description.y")
    POI_phosPeps<- POI_matrix[,1:(idx_split-1)]
    POI_globProts<- POI_matrix[,(idx_split-1):ncol(POI_matrix)]
    # Do the plotting
    .plotProteinNPhosphoPeptidesNtoC_FragPipe_TMT(protName = MyProteinName_poi, globalProtFC = POI_protFC,
                                              POI_tableProtein = POI_globProts, POI_tablePeptide = POI_phosPeps)
  }
  dev.off()
}




#' write and render pdfs for phospho Rmarkdown
#' @param grp2 grp2 object with all contrasts and all results
#' @param name experiment name refound in zip and htmls
#' @param ZIPDIR directory name to write results in
#' @param boxplot boolean to either do all boxplots or not
#' @export
#'
write_phosphoDEA_all <- function (grp2, name, ZIPDIR, boxplot = TRUE)
{
  fname <- paste0("DE_", name)
  qcname <- paste0("QC_", name)
  outpath <- file.path(ZIPDIR, fname)
  logger::log_info("writing into : ", outpath, " <<<<")
  # write result files ORA, GSEA, xlsx
  prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname)

  prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname, markdown = "_Grp2Analysis_Phospho.Rmd")
  prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname,
                         markdown = "_DiffExpQC_Phospho.Rmd")
  bb <- grp2$RES$transformedlfqData
  grsizes <- dplyr::pull(dplyr::summarize(dplyr::group_by(bb$factors(),
                                                          dplyr::across(bb$config$table$factor_keys_depth())),
                                          n = n()), n)
  if (boxplot) {
    if (sum(!grepl("^control", bb$config$table$factor_keys(),
                   ignore.case = TRUE)) > 1 & all(grsizes == 1)) {
      prolfquapp::writeLinesPaired(bb, outpath)
    }
    else {
      pl <- bb$get_Plotter()
      pl$write_boxplots(outpath)
    }
  }
}

# MS-stats like normalization for protein change
# https://bioc.ism.ac.jp/packages/3.12/bioc/vignettes/MSstatsTMTPTM/inst/doc/MSstatsTMTPTM.Workflow.html
# look for: how to adjust PTMs

#' Apply MSstatsPTM like site normalization (adjustment) for the protein fold-change and pvalues
#' @param mycombo dataframe from PTMsite and protein prolfqua statistics
#' @export
#'
doMSstatsLikeSiteNormalizationUsingProteinStatsOnComboObject <- function (mycombo)
{
  resultCombo <- data.frame(stringsAsFactors = TRUE)
  for (i in 1:length(unique(mycombo$contrast))) {
    OneC <- mycombo[mycombo$contrast == unique(mycombo$contrast)[i],]
    OneC$MSstatsPTMadj_log2fc <- OneC$diff.x - OneC$diff.y
    OneC$MSstatsPTMadj_s2 <- OneC$std.error.x^2
    OneC$MSstatsPTMadj_s2prot <- OneC$std.error.y^2
    OneC$MSstatsPTMadj_stderr <- sqrt(OneC$MSstatsPTMadj_s2 + OneC$MSstatsPTMadj_s2prot)
    OneC$MSstatsPTMadj_numer <- (OneC$MSstatsPTMadj_s2 + OneC$MSstatsPTMadj_s2prot)^2
    OneC$MSstatsPTMadj_denom <- (OneC$MSstatsPTMadj_s2^2 / OneC$df.x + OneC$MSstatsPTMadj_s2prot^2 / OneC$df.y)
    OneC$MSstatsPTMadj_df <- OneC$MSstatsPTMadj_numer / OneC$MSstatsPTMadj_denom
    OneC$MSstatsPTMadj_tval <- OneC$MSstatsPTMadj_log2fc / OneC$MSstatsPTMadj_stderr
    OneC$MSstatsPTMadj_pVals <- 2 * stats::pt(abs(OneC$MSstatsPTMadj_tval), OneC$MSstatsPTMadj_df, lower.tail = FALSE)
    #adjust pV for multiple testing
    OneC$MSstatsPTMadj_FDR <- p.adjust(OneC$MSstatsPTMadj_pVals)
    resultCombo <- rbind(resultCombo, OneC)
  }
  return(resultCombo)
}


# ////////// to work on

#' Generate differential expression analysis reports
#'
#' Writes results of DEA see \code{\link{generate_DEA_reports}}
#' @export
#'
write_DEA_all <- function(grp2, name, ZIPDIR, boxplot = TRUE){
  fname <- paste0("DE_", name)
  qcname <- paste0("QC_", name)
  outpath <- file.path( ZIPDIR, fname)
  logger::log_info("writing into : ", outpath, " <<<<")
  prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname)
  prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname)
  prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname, markdown = "_DiffExpQC.Rmd")

  bb <- grp2$RES$transformedlfqData
  grsizes <- bb$factors() |>
    dplyr::group_by(dplyr::across(bb$config$table$factor_keys_depth())) |>
    dplyr::summarize(n = n()) |>
    dplyr::pull(n)
  if (boxplot) {
    if (sum(!grepl("^control",bb$config$table$factor_keys(), ignore.case = TRUE))  > 1 &
        all(grsizes == 1)
    ) {
      prolfquapp::writeLinesPaired(bb, outpath)
    } else {
      pl <- bb$get_Plotter()
      pl$write_boxplots(outpath)
    }
  }
}



#' Write differential expression analysis results
#'
#' @rdname make_DEA_report
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @param outpath path to place output
#' @param xlsxname file name for xlsx
#' @export
#' @family workflow
#'
write_DEA <- function(GRP2, outpath, xlsxname = "AnalysisResults"){
  dir.create(outpath)
  rd <- GRP2$RES$lfqData
  tr <- GRP2$RES$transformedlfqData
  ra <- GRP2$RES$rowAnnot
  formula <- data.frame(
    formula = GRP2$RES$formula,
    contrast_name = names(GRP2$pop$Contrasts),
    contrast = GRP2$pop$Contrasts)

  wideraw <- dplyr::inner_join(ra$row_annot, rd$to_wide()$data, multiple = "all")
  widetr <- dplyr::inner_join(ra$row_annot , tr$to_wide()$data, multiple = "all")

  ctr <- dplyr::inner_join(ra$row_annot , GRP2$RES$contrMerged$get_contrasts(), multiple = "all")
  resultList <- list()
  resultList$annotation = tr$to_wide()$annot
  resultList$normalized_abundances = dplyr::inner_join(ra$row_annot, tr$data,multiple = "all")
  resultList$raw_abundances_matrix = wideraw
  resultList$normalized_abundances_matrix = widetr
  resultList$diff_exp_analysis = ctr
  resultList$formula = formula
  resultList$summary = GRP2$RES$Summary
  resultList$missing_information = prolfqua::UpSet_interaction_missing_stats(rd$data, rd$config, tr = 1)$data

  # add protein statistics
  st <- GRP2$RES$transformedlfqData$get_Stats()
  resultList$protein_variances <- st$stats()

  bkg <- GRP2$RES$rowAnnot$row_annot$IDcolumn
  ff <- file.path(outpath ,"ORA_background.txt")
  write.table(bkg,file = ff, col.names = FALSE,
              row.names = FALSE, quote = FALSE)

  fg <- GRP2$RES$contrastsData_signif
  ora_sig <- split(fg$IDcolumn, fg$contrast)

  for (i in names(ora_sig)) {
    ff <- file.path(outpath, paste0("Ora_",i,".txt" ))
    logger::log_info("Writing File ", ff)
    write.table(ora_sig[[i]],file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE)
  }

  fg <- GRP2$RES$contrastsData
  gsea <- fg |> dplyr::select( contrast, IDcolumn, statistic) |> dplyr::arrange( statistic )
  gsea <- split(dplyr::select( gsea, IDcolumn, statistic ), gsea$contrast)

  for (i in names(gsea)) {
    ff <- file.path(outpath, paste0("GSEA_",i,".rnk" ))
    logger::log_info("Writing File ", ff)
    write.table(na.omit(gsea[[i]]),file = ff, col.names = FALSE,
                row.names = FALSE, quote = FALSE, sep = "\t")
  }
  if (nrow(resultList$normalized_abundances) > 1048575) {
    resultList$normalized_abundances <- NULL
  }
  writexl::write_xlsx(resultList, path = file.path(outpath, paste0(xlsxname, ".xlsx")))
}

#' Render DEA analysis report
#' @rdname make_DEA_report
#' @param GRP2 return value of \code{\link{make_DEA_report}}
#' @param outpath path to place output
#' @param htmlname name for html file
#' @param word default FALSE, if true create word document.s
#' @param markdown which file to render
#' @export
#' @family workflow
render_DEA <- function(GRP2,
                       outpath,
                       htmlname="Result2Grp",
                       word = FALSE,
                       markdown = "_Grp2Analysis.Rmd"){
  dir.create(outpath)

  rmarkdown::render(
    markdown,
    params = list(grp = GRP2) ,
    output_format = if (word) {
      bookdown::word_document2(toc = TRUE, toc_float = TRUE) } else {
        bookdown::html_document2(toc = TRUE, toc_float = TRUE)
      }
  )
  fname <- paste0(tools::file_path_sans_ext(markdown), if (word) {".docx"} else {".html"})
  if (file.copy(fname, file.path(outpath, paste0(htmlname,if (word) {".docx"} else {".html"})), overwrite = TRUE)) {
    file.remove(fname)
  }
}


# Do the NtoCplot
# 2023-07-05 adapted for FP outputs
# this undocumented function is used to plot the phospho peptides from NtoC-term with indicated log2-fold-change-bars
# the function is used inside another documented function -> the 2019 function should be used from 2019 onwards
.plotProteinNPhosphoPeptidesNtoC_FragPipe_TMT <- function(protName, globalProtFC, POI_tableProtein, POI_tablePeptide,
                                                      protColor = "blue", protWtdh = 4, maxX = 1200, protNCoffset = 40,
                                                      pepSigStarOffset = 0.05, pepSigStarSize = 2, fdrRelaxedThreshold = 0.2,
                                                      fdrStringentThreshold = 0.05, greyLineLog2Threshold = 1) {
  # issue: we do have some not localized sites in, filter these but keep original
  POI_pepOriginal <- POI_tablePeptide
  POI_tablePeptide <- POI_tablePeptide[POI_tablePeptide$SinglePhosLocalized_bool,]

  # handle cases wher we only have non localized peptides

  myYaxisLimiter <- max(ceiling(max(abs(POI_tablePeptide$diff.x))), ceiling(max(abs(POI_tableProtein$diff.y))), na.rm = TRUE)
  miny <- -myYaxisLimiter
  maxy <- myYaxisLimiter
  # get positions
  Phos_PositionsInProteins <- vector(length=nrow(POI_tablePeptide[POI_tablePeptide$SinglePhosLocalized_bool,]))
  for (i in 1:nrow(POI_tablePeptide)) {
    Phos_PositionsInProteins[i] <-  as.numeric(POI_tablePeptide$posInProtein[i])
  }
  # careful here we do not want positions that are parsed from multiply phos
  scaleX <- 1/max(Phos_PositionsInProteins, na.rm = TRUE)*1200

  # handle if no peptide is localized
  if (nrow(POI_tablePeptide) != 0)  POI_tablePeptide$PositionsINproteins <- Phos_PositionsInProteins

  phosphoPlotTitle <- paste("Prot: ",protName,"\n Peptides from: ", POI_pepOriginal$protein_Id[1],"# phospho = ",
                            length(POI_pepOriginal$diff.x), sep=" ")

  # Do the plot start with white
  plot(c(-100,0,1400, 0,0), c(globalProtFC,0,0,miny, maxy), pch=".", col="white", main=phosphoPlotTitle, ylab="log2FC",
       xlab="Full length protein", xaxt='n')
  axis(side = 1,at = c(-100,0,1200),labels = c("protFC","0","100%"))
  #grey box for prot FC
  rect(xleft = -150,ybottom = -100,xright = -50,ytop = 100, col = "lightgrey",border = FALSE)
  if (is.na(globalProtFC)) {
    points(-100,0,pch="0", col="orange", cex=3)
  } else {
    segments(-100,0,-100, globalProtFC, col=protColor, lwd=protWtdh)
  }

  # set the stage for the protein
  segments(0,0,maxX,0, col="black", lwd=1.5)
  points(-protNCoffset,0,pch="N")
  points(maxX+protNCoffset,0, pch="C")

  # plot all peptides if we have something in peptide table
  if (nrow(POI_tablePeptide) != 0) {
    #make modAA as factor
    POI_tablePeptide$modAA <- as.factor(POI_tablePeptide$AA)

    #make modelName a factor for lty -> strange error
    POI_tablePeptide$modelName.x <- as.factor(POI_tablePeptide$modelName.x)

    for (i in 1:length(POI_tablePeptide$PositionsINproteins)) {
      if (!is.na(POI_tablePeptide$PositionsINproteins[i])) {
        # plot lines
        if(unique(POI_tablePeptide$modelName.x) == "Linear_Model_moderated") {
          segments(POI_tablePeptide$PositionsINproteins[i]*scaleX, 0,
                   POI_tablePeptide$PositionsINproteins[i]*scaleX, POI_tablePeptide$diff.x[i],
                   col=as.numeric(POI_tablePeptide$modAA[i]), lwd=1, lty=1) # hack: when all same model = lty = 1
        }
        else {
          segments(POI_tablePeptide$PositionsINproteins[i]*scaleX, 0,
                   POI_tablePeptide$PositionsINproteins[i]*scaleX, POI_tablePeptide$diff.x[i],
                   col=as.numeric(POI_tablePeptide$modAA[i]), lwd=1, lty=abs(as.numeric(POI_tablePeptide$modelName.x[i])-2)+1) # quite a hack for lty to get 2 = 1
        }
      }


      # decoration!
      # site is stringent significant
      if (POI_tablePeptide$FDR.x[i] < fdrStringentThreshold && !is.na(POI_tablePeptide$modAA[i])) {
        if(POI_tablePeptide$diff.x[i] > 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                    scaleX, POI_tablePeptide$diff.x[i]+pepSigStarOffset, pch="*",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(POI_tablePeptide$modAA[i]))

        if(POI_tablePeptide$diff.x[i] < 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                    scaleX, POI_tablePeptide$diff.x[i]-pepSigStarOffset, pch="*",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(POI_tablePeptide$modAA[i]))
      }
      else if (POI_tablePeptide$FDR.x[i] < fdrRelaxedThreshold && !is.na(POI_tablePeptide$modAA[i])) {
        if(POI_tablePeptide$diff.x[i] > 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                    scaleX, POI_tablePeptide$diff.x[i]+pepSigStarOffset, pch="+",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(POI_tablePeptide$modAA[i]))

        if(POI_tablePeptide$diff.x[i] < 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                    scaleX, POI_tablePeptide$diff.x[i]-pepSigStarOffset, pch="+",
                                                  cex=pepSigStarSize,
                                                  col=as.numeric(POI_tablePeptide$modAA[i]))
      }
      else if (POI_tablePeptide$FDR.x[i] > fdrRelaxedThreshold && !is.na(POI_tablePeptide$modAA[i])) {
        if(POI_tablePeptide$diff.x[i] > 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                    scaleX, POI_tablePeptide$diff.x[i]+pepSigStarOffset, pch="x",
                                                  cex=pepSigStarSize/3,
                                                  col=as.numeric(POI_tablePeptide$modAA[i]))

        if(POI_tablePeptide$diff.x[i] < 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                    scaleX, POI_tablePeptide$diff.x[i]-pepSigStarOffset, pch="x",
                                                  cex=pepSigStarSize/3,
                                                  col=as.numeric(POI_tablePeptide$modAA[i]))
        if(POI_tablePeptide$diff.x[i] == 0) points(POI_tablePeptide$PositionsINproteins[i] *
                                                     scaleX, POI_tablePeptide$diff.x[i], pch="+",
                                                   cex=pepSigStarSize/3,
                                                   col=as.numeric(POI_tablePeptide$modAA[i]))
      }

    }
  }

  # add legend for not localized phospho
  legend("bottomleft", legend = paste("Not fully localized phospho peptides:", sum(POI_pepOriginal$SinglePhosLocalized_bool == FALSE)))
  # no phospho site is localized
  if (length(unique(na.omit(POI_tablePeptide$modAA))) != 0) {
    legend("bottomright", legend = unique(na.omit(POI_tablePeptide$modAA)),
           col=unique(as.numeric(na.omit(POI_tablePeptide$modAA))),
           text.col = unique(as.numeric(POI_tablePeptide$modAA)), lwd=1)
    legend("topright", legend = c("dashed line -> imputed model", "below stringent FDR threshold", "below relaxed FDR threshold", "above FDR threshold"),
           col = c("orange", "black", "black", "black"), pch=c(" ","*", "+", "x"), cex = 0.7, text.col = c("orange", "black", "black", "black"))
  }

  abline(h=greyLineLog2Threshold, col="grey", lty=2)
  abline(h=-greyLineLog2Threshold, col="grey", lty=2)
}

