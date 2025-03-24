#' read AA fasta file, identify decoy entries and return fw fasta object only
#' @param FastaFileName filename of fasta file
#' @param decoyPattern pattern for decoy accessions
#' @return fastaObject without decoy sequences
#' @export
#'
read_fasta <- function(FastaFileName, decoy_pattern = "^REV_") {
  myFasta_decoy <- seqinr::read.fasta(file = FastaFileName, seqtype = "AA", as.string = TRUE)
  seqNames <- seqinr::getName(myFasta_decoy)
  nodeoySeq <- myFasta_decoy[!grepl( decoy_pattern , seqNames )]
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
  if ((soi - seqWindowOffset) < 0) {
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
    #rm(POI_matrix)
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



# For multiple plex TMT experiments each plex writes a separate psm file
# They shall be read in individually till long format (using tidy_FragPipe_psm)
# and then combined and  passed to this function here
#

# To do: write docu to get this function exported
preprocess_FP_multiplexPSM <- function (psm, fasta_file, annotation,
                                        purity_threshold = 0.5, PeptideProphetProb = 0.9,
                                        column_before_quants = c("Quan Usage","Mapped Proteins"),
                                        pattern_contaminants = "^zz|^CON",  pattern_decoys = "rev_"){
  annot <- annotation$annot
  atable <- annotation$atable
  annot <- dplyr::mutate(annot, raw.file = gsub("^x|.d.zip$|.raw$",
                                                "", (basename(annot[[atable$fileName]]))))
  #psm <- prolfquapp::tidy_FragPipe_psm(quant_data, column_before_quants = column_before_quants) # we do this step outside
  nrPeptides_exp <- psm$nrPeptides
  psm <- psm$data
  psm$qValue <- 1 - psm$Probability
  nr <- sum(annot[[annotation$atable$fileName]] %in% sort(unique(psm$channel)))
  logger::log_info("nr : ", nr, " files annotated out of ",
                   length(unique(psm$channel)))
  stopifnot(nr > 0)

  logger::log_info("channels in annotation which are not in psm.tsv file : ",
                   paste(setdiff(annot[[annotation$atable$fileName]], sort(unique(psm$channel))),
                         collapse = " ; "))
  logger::log_info("channels in psm.tsv which are not in annotation file : ",
                   paste(setdiff(sort(unique(psm$channel)), annot[[annotation$atable$fileName]]),
                         collapse = " ; "))
  atable$ident_Score = "Probability"
  atable$ident_qValue = "qValue"
  atable$hierarchy[["protein_Id"]] <- c("Protein")
  atable$hierarchy[["peptide_Id"]] <- c("Peptide")
  atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide",
                                            "Assigned.Modifications")
  atable$set_response("abundance")
  bycol <- c("channel")
  names(bycol) <- atable$fileName
  psma <- dplyr::inner_join(annot, psm, multiple = "all", by = bycol)
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(psma, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
  fasta_annot <- dplyr::left_join(nrPeptides_exp, fasta_annot,
                                  by = c(Protein = "fasta.id"))
  fasta_annot <- dplyr::rename(fasta_annot, `:=`(!!lfqdata$config$table$hierarchy_keys_depth()[1],
                                                 !!sym("Protein")))
  fasta_annot <- dplyr::rename(fasta_annot, description = fasta.header)
  prot_annot <- prolfquapp::ProteinAnnotation$new(lfqdata,
                                                  fasta_annot, description = "description", cleaned_ids = "proteinname",
                                                  full_id = "protein_Id", exp_nr_children = "nrPeptides",
                                                  pattern_contaminants = pattern_contaminants, pattern_decoys = pattern_decoys)
  lfqdata$remove_small_intensities()
  return(list(lfqdata = lfqdata, protein_annotation = prot_annot))
}


