# Author: Lauren Rylaarsdam, PhD
# 2024-2025

############################################################################################################################
#' @title amethyst-class
#' @description An S4 class to store and manipulate single-cell methylation data
#'
#' @slot h5paths Path to the hdf5 file containing base-level read information organized by methylation type and barcode.
#' If using Facet, this file should also contain aggregated methylation observations over larger genomic features
#' @slot genomeMatrices Slot to store aggregated genomic information
#' @slot tracks Aggregated methylation levels over short genomic regions. For plotting purposes (e.g. histograM, heatMap)
#' @slot reductions Slot to store dimensionality reductions, such as irlba and umap/tsne
#' @slot index Corresponding chromosome coordinates for each cell in the hdf5 file
#' @slot metadata Optional cell metadata. If included, make sure row names are cell IDs
#' @slot ref Genome annotation file with chromosome, start, and end position information for genes of interest. See the "makeRef" function
#' @slot results Convenience slot to store results if desired
#'
#' @import methods
#' @export
methods::setClass("amethyst", slots = c(
  h5paths = "ANY",
  genomeMatrices = "ANY",
  tracks = "ANY",
  reductions = "ANY",
  index = "ANY",
  metadata = "ANY",
  ref = "ANY",
  results = "ANY"
))

############################################################################################################################
#' @title createObject
#' @description Create object of class amethyst
#'
#' @param h5paths Path to the hdf5 file containing base-level read information organized by methylation type and barcode.
#' If using Facet, this file should also contain aggregated methylation observations over larger genomic features
#' @param index Corresponding chromosome coordinates for each cell in the hdf5 file
#' @param metadata Optional cell metadata. If included, make sure row names are cell IDs
#' @param genomeMatrices Slot to store aggregated genomic information
#' @param reductions Slot to store dimensionality reductions, such as irlba and umap/tsne
#' @param ref Genome annotation file with chromosome, start, and end position information for genes of interest. See the "makeRef" function
#' @param tracks Aggregated methylation levels over short genomic regions. For plotting purposes (e.g. histograM, heatMap)
#' @param results Convenience slot to store results if desired
#'
#' @return Returns a single object of class amethyst with h5path, genomeMatrices, reduction, index, metadata, and ref slots.
#' @importFrom methods new
#' @export
#' @examples obj <- createObject()

createObject <- function(h5paths = NULL,
                         genomeMatrices = NULL,
                         tracks = NULL,
                         reductions = NULL,
                         index = NULL,
                         metadata = NULL,
                         ref = NULL,
                         results = NULL) {
  methods::new(Class = "amethyst",
               h5paths = h5paths,
               genomeMatrices = genomeMatrices,
               tracks = tracks,
               reductions = reductions,
               index = index,
               metadata = metadata,
               ref = ref,
               results = results)
}

############################################################################################################################
# With input from Dave Ross and Felix Schlesinger at Scale Biosciences
#' @title createScaleObject
#' @description Helper function for converting Scale Biosciences pipeline output into an Amethyst object
#'
#' @param directory Path to the directory containing sample output from the Scale Biosciences computational pipeline.
#' Path should be structured as follows: "path_to_folder/ScaleMethyl.out/samples". The folder is expected to contain
#' an ".allCells.csv" file with metadata, a genome_bin_matrix folder with pre-constructed matrices, and a
#' methylation_coverage folder containing .h5 files with base-level methylation information for each cell.
#' @param genomeMatrices Optional name of pre-constructed matrices in the genome_bin_matrix folder to include.
#' @return Returns a populated amethyst object for futher analysis.
#' @importFrom data.table fread setnames rbindlist
#' @importFrom dplyr filter bind_rows
#' @importFrom utils read.csv
#' @importFrom Matrix readMM
#' @export
#' @examples
#' \dontrun{
#'   obj <- createScaleObject(directory = "~/Downloads/ScaleMethyl.out/samples", genomeMatrices = list("CG.score", "CH"))
#' }
createScaleObject <- function(directory,
                              genomeMatrices = NULL) {
  samples <- sub("\\.allCells\\.csv$", "", grep("\\.allCells\\.csv$", list.files(directory), value = TRUE))

  metadata <- list()
  for (i in samples) {
    metadata[[i]] <- read.csv(paste0(directory, "/", i, ".allCells.csv")) |> dplyr::filter(pass == "pass")
  }

  h5paths <- list()
  for (i in samples) {
    h5paths[[i]] <- data.frame(barcode = metadata[[i]]$cell_id,
                               path = file.path(directory, "methylation_coverage", "amethyst", i, paste0(i, ".", metadata[[i]]$tgmt_well, "_cov.h5")))
  }

  obj <- createObject(
    metadata = data.table::rbindlist(metadata) |> tibble::column_to_rownames(var = "cell_id"),
    h5paths = dplyr::bind_rows(h5paths),
  )

  if (!is.null(genomeMatrices)) {
    mtx <- list()
    for (i in genomeMatrices) {
      for (j in samples) {
        mtx_path <- paste0(directory, "/genome_bin_matrix/", j, ".", i, ".mtx.gz")
        if(length(k <- grep("CG", i))) {
          feature_path <- paste0(directory, "/genome_bin_matrix/", j, ".CG.features.tsv")
        } else if (length(k <- grep("CH", i))) {
          feature_path <- paste0(directory, "/genome_bin_matrix/", j, ".CH.features.tsv")
        }
        features <- data.table::fread(feature_path, header = FALSE)
        data.table::setnames(features, "features")

        mtx[[i]][[j]] <- as.array(Matrix::readMM(mtx_path))
        dimnames(mtx[[i]][[j]]) <- list(features$features, metadata[[j]]$cell_id)
        mtx[[i]][[j]] <- as.data.frame(mtx[[i]][[j]]) |> tibble::rownames_to_column(var = "window")
      }
      mtx[[i]] <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), mtx[[i]]) |> tibble::column_to_rownames(var = "window")
    }
    obj@genomeMatrices <- mtx
  }
  return(obj)
}

############################################################################################################################
#' @title extractAttributes
#' @description Extract information from the attributes column of a gtf file
#' This helper function was taken from https://www.biostars.org/p/272889/
#'
#' @param gtf_attributes A string containing the attributes column from a GTF file
#' @param att_of_interest Name of the attribute containing the desired information
#' @return Returns the value associated with the specified attribute key from the GTF attributes string
#' @export
#' @examples
#' \dontrun{
#'   gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))
#' }
extractAttributes <- function(gtf_attributes,
                              att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return(unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

############################################################################################################################
#' @title makeRef
#' @description Helper function to generate a genome annotation file in the expected format
#'
#' @param genome Reference genome, e.g. "hg38"
#' @param gtf If not using hg19, hg38, mm10, or mm39, please provide the gtf.gz file path
#' @param attributes Information to extract from the gtf file. Must be a column name
#' @return Returns an annotated reference of gene locations
#' @importFrom dplyr mutate
#' @importFrom rtracklayer readGFF
#' @export
#' @examples
#' \dontrun{
#'   obj@ref <- makeRef(genome = "hg38")
#' }
makeRef <- function(genome,
                    gtf = NULL,
                    attributes = c("gene_name", "exon_number")) {
  if (genome == "hg38") {
    options(timeout = 1000)
    gtf <- rtracklayer::readGFF("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz")
  } else if (genome == "hg19") {
    options(timeout = 1000)
    gtf <- rtracklayer::readGFF("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz")
  } else if (genome == "mm10") {
    options(timeout = 1000)
    gtf <- rtracklayer::readGFF("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz")
  } else if (genome == "mm39") {
    options(timeout = 1000)
    gtf <- rtracklayer::readGFF("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz")
  } else if (!is.null(gtf)) {
    gtf <- rtracklayer::readGFF(gtf)
  } else {
    stop("Ref was different from defaults but not provided. Please provide unzipped gtf reference annotation file.")
  }
  for (i in attributes) {
    gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))
  }
  gtf <- gtf |> dplyr::mutate(location = paste0(seqid, "_", start, "_", end))
  gtf <- data.table::data.table(gtf)

  return(gtf)
}

############################################################################################################################
#' @title fetchMarkers
#' @description Fetch a pre-selected list of key marker genes
#' Key marker genes for brain methylation selected from http://neomorph.salk.edu/omb/ct_table
#'
#' @param species Human or mouse. Determines gene case.
#' @param type Tissue type of interest. For now, "brain" and "pbmc" are options; others will be added
#' @return Returns a character vector of key marker genes.
#' @export
#' @examples
#' \dontrun{
#'   markerGenes <- fetchMarkers(species = "human", type = "brain")
#' }
fetchMarkers <- function(species,
                         type) {
  if (type == "brain") {
    markerGenes <- c("Abat", "Abca12","Abca13","Abhd2","Ablim1","Ablim2", "Ablim3","Acan","Acsf3","Adamts15", "Adcy8","Adcyap1r1","Adgrf5","Adgrg4", "Adgrl2","Adora2a","Adra1a","Afap1","Aff2", "Ajap1","Ak5","Akap13","Ankk1","Ano1", "Aplnr","Aqp4","Arhgap15","Arhgap22","Arhgap25",
                     "Arhgef28","Arid5b","Arl15","Asap1","Asb2", "Astn2","Atp1a2","Atp2b4","Auts2","B3gat1", "B3gat2","Babam2","Bcl11b", "Bcl11a","Bcl2l11", "Bend5","Bend7","Bmp7","Bmper","Bmpr1b", "Bnc2","Boc","Brinp3","Brsk2","C2cd3", "Cachd1","Cacna1e","Cacna1h","Cacna1i",
                     "Cacna2d1","Cacna2d2","Cacng3","Cacng4", "Cacng5","Cadm1","Cadm2","Cadps2","Calb2", "Cald1","Caln1","Camta1","Car2","Carmil1","Cars","Casz1","Cbfb","Ccbe1","Ccdc88c", "Ccser1","Cd9","Cdh13","Cdh18","Cdh2", "Cdh22","Cdh4","Cdh6","Cdh7","Cdh8","Cdh9",
                     "Cdk14","Cdk17","Cdk19","Cdyl2","Celf4", "Cemip","Cflar","Chat","Chchd6","Chd7", "Chodl","Chrm3","Chrna7","Chst11","Clic4", "Clic5","Clmn","Clmp","Clstn2","Cnih3", "Cntn4","Cntn6","Cntnap4","Cntnap5a", "Cntnap5b","Cntnap5c","Col12a1","Col14a1",
                     "Col15a1","Col4a1","Col6a3","Colgalt2","Cpa6", "Cpe","Cpeb4","Crim1","Crot","Cryz12", "Csf1r","Csf2rb2","Csgalnact1","Csmd1", "Csmd2","Ctnna2","Cux1","Cux2","Cyp7b1", "Daam2","Dab1","Dchs2","Ddah1","Dennd1a", "Deptor","Dgki","Dlc1","Dlg5","Dlgap1",
                     "Dnah9","Dnajc15","Dner","Dock10","Dock4", "Dock5","Dock9","Dok5","Dpf3","Dpy19l1", "Drd1","Drd2","Dscam","Dscaml1","Dusp10","Ecel1","Efnb2","Egf17","Egflam","Egfr", "Egln3","Elavl2","Elavl4","Elfn2", "Enox1","Enpp2","Enpp6","Entpd3","Epha3",
                     "Epha4","Epha5","Epha6","Epha7","Ephb2", "Eps8","Epsti1","Erg","Eri3","Erich3","Ermp1", "Ern1","Esrrg","Etv1","Etv6","Evc2","Eya1","Eya2","Fabp7","Fam107a","Fam129a","Fam174a", "Fam196b","Fam222a","Fat1","Fat3", "Fat4","Fbn1","Fcrls","Fgd5","Fgf13","Fibcd1",
                     "Flrt2","Fmnl2","Fn1","Foxn3","Foxo3", "Foxp1","Foxp2","Fras1","Frmd3","Frmd6", "Fstl4","Fyn","Fzd3","Gabbr2","Gabra2","Gad1","Gad2","Galnt10","Galnt17","Galnt9", "Galntl6","Gap43","Garnl3","Gas7","Gcn1l1","Gfap","Gfra1","Gfra2","Glant17","Gm10421",
                     "Gm11267","Gm28308","Gnal","Gng12", "Golim4","Gpd2","Gpr26","Grid1","Grik3", "Grik4","Grin3a","Grk5","Grm1","Grm3","Grm4","Grm7","Grm8","Gse1","Gucy1b1","H2afy2", "Hcrtr2","Hdac4","Heg1","Hipk2","Homer2", "Hpse","Hrh1","Hs3st4","Htr2a","Htr7",
                     "Hunk","Igdcc3","Igf1r","Igsf3","Ikzf1", "Il1rap","Il1rl2","Ildr2","Inpp4b", "Insyn2a","Iqcm","Itgb8","Jup","Kank1","Kcnab1","Kcnc2","Kcnd3","Kcnh1","Kcnh5", "Kcnip2","Kcnip4","Kcnj6","Kcnk10","Kcnk2", "Kcnk3","Kcnma1","Kcnmb4","Kcnn3","Kctd8",
                     "Kdr","Khdrbs2","Khdrbs3","Kif13a","Kif26a", "Kif5c","Kit","Kitl","Klh12","Ksr1", "L3mbtl4","Lama4","Lamb1","Lamp5","Laptm5","Lats2","Ldb2","Ldlrad4","Lhfp","Lhx2", "Lhx6","Lima1","Limch1","Lingo1","Lingo2","Lmo1","Lnx2","Lpp","Lrfn2","Lrp1b",
                     "Lrrc4","Lrrk1","Lrrtm3","Lrrtm4","Lypd1", "Lypd6","Lypd6b","Magi1","Man1a","Man1c1","Man2a1","Map3k1","Map3k21","Map3k5","Map4","Map4k4","Mapk10","Mapk4","March4", "Masp1","Matn2","Mbp","Mcc","Mctp2", "Med12l","Mef2c","Megf11","Megf9","Meis1",
                     "Meis2","Mertk","Metrn","Mgat5","Mgll","Mob2", "Mobp","Mog","Morn5","Mpdz","Mpped1", "Mthfd1l","Mtus2","Myh13","Myo16","Myo1b","Myocd","Ncam2","Ndrg1","Ndst3","Ndst4", "Nectin1","Nectin3","Nell1","Neurl1a","Neurod1","Neurod2","Neurod6","Nfasc","Nfia","Nfib","Nfix","Nkd1",
                     "Nkd2","Nos1","Nos1ap","Notch2","Npas2", "Npas3","Npsr1","Nr3c2","Nr4a2","Nr4a3", "Nrip","Nrp1","Nrp2","Nrxn1","Nrxn3", "Ntf3","Ntn1","Ntng1","Ntng2","Nuak1","Nxn", "Nxph1","Ofcc1","Olig1","Olfml2b","Onecut2", "Opn3","Oprd1","Osbpl3","Oxr1","P2ry13",
                     "P3h2","Pacs1", "Pag1","Pakap","Pam","Pappa","Pappa2", "Paqr8","Pard3","Pard3b","Pard6g", "Parm1","Parp8","Pax6","Pbx1","Pbx3","Pcdh19","Pcdh7","Pcsk5","Pcsk6","Pde4b","Pde4d", "Pde8b","Pdgfra","Pdgfrb","Pdia5","Peak1", "Pear1","Pexl4","Pfkfb3","Phactr4","Phf14",
                     "Pip4k2a","Pip5k1b","Pkb4","Pkhd1", "Pknox2","Plcb4","Plch2","Plcl1","Plcxd3", "Pld5","Plekha2","Plekha6","Plekhg1","Plpp3","Plpp4","Plppr1","Plxdc1","Plxdc2", "Plxna1","Plxna4","Plxnc1","Pmp22","Pou6f2", "Ppargc1a","Ppfibp1","Ppm1l","Ppp1r12b",
                     "Ppp1r14c","Prdm16","Prex1","Prex2","Prima1", "Prkag2","Prkcb","Prkch","Prkcq", "Prkd1","Prlr","Prox1","Prr5l","Ptbp3","Ptgfrn","Ptn","Ptpn3","Ptpn4","Ptprd","Ptpre", "Ptprf","Ptprg","Ptprk","Ptprm","Ptpro", "Ptprr","Ptprt","Ptpru","Ptprz1","Pvalb",
                     "Qk","Rab11fip2","Rab3b","Rab3c", "Ralgapa2","Rapgef5","Rarb","Rasgef1b","Rasgrf2","Rbfox3","Rbms1","Rdh10","Reps2","Rerb1","Rerg","Rffl","Rfx3","Rgl1","Rgs12","Rgs5", "Rgs8","Ripor2","Rlen","Rmst","Rnf144a", "Rnf144b","Rnf152","Robo1","Robo2","Rora",
                     "Rorb","Rph3a","Rspo2","Runx1t1","Rxra","S100b","Sall3","Satb1","Satb2","Scube1","Sel1l3", "Sema3e","Sema3g","Sema4d","Sema5a","Sema5b", "Sema6a","Sema6d","Sema7a","Serpine2", "Sesn3","Sgcd","Sgk1","Sh3pxd2a","Sh3pxd2b", "Shb","Shc3","Sik3","Sipa1l2","Slc17a6","Slc17a7","Slc17a8",
                     "Slc1a2","Slc1a3","Slc18a3","Slc24a2","Slc24a4", "Slc35f4","Slc39a8","Slc44a5","Slc4a4", "Slc5a7","Slc6a1","Slc6a13","Slc7a5","Slc8a1","Slc8a3","Slit1","Slit3","Smad4","Smad6", "Smad7","Sntb1","Sntg2","Snx7","Sobp", "Sorcs1","Sorcs2","Sox13","Sox2","Sox6","Spats2l",
                     "Specc1","Spock3","Spon1","Srgap1", "Ssbp3","Sst","St18","St3gal1","St3gal6", "St5","St6galnac3","St6galnac5","Stat5b","Stk32b","Stk39","Strip2","Stxbp6","Sulf2", "Susd1","Susd4","Susd5","Sv2b","Svil", "Syn3","Syndig1","Synpr","Syt6","Syt7",
                     "Tbc1d1","Tbc1d4","Tbr1", "Tcerg1l","Tcf4","Tcf7l2", "Tead1","Tenm2","Tenm3","Tenm4","Tex9", "Tgfa","Tgm3","Th","Thsd7a","Tle4","Tll1", "Tmem132b","Tmem132c","Tmem132d","Tmem132e", "Tmem178","Tmem178b","Tmem65","Tmtc1", "Tmtc2","Tnks","Tnr","Tns3","Tox","Tox3",
                     "Trabd2b","Trim2","Trim9","Trpc3","Trpc4", "Trpc7","Trps1","Tshz1","Tshz2","Tshz3", "Tspan11","Tspan14","Tspan5","Tspan9","Ttc39c","Ttr","Tunar","Ubn2","Ubtd1","Unc13c", "Unc5b","Unc5c","Uox","Upp1","Usp6nl", "Vat1l","Vav2","Vgll4","Vip","Vps13d",
                     "Vstm2b","Vwc2l","Whamm","Whrn","Wnt3", "Wnt9b","Wwc2","Wwp2","Xkr6","Xkr7","Xpr1", "Xxylt1","Xylt1","Zdhhc14","Zeb2","Zfand4","Zfp366","Zfp385b","Zfp423","Zfp462", "Zfp536","Zfp618","Zfp710","Zfp804a", "Zfp827","Zfpm2","Zhx2","Zmat4","Zmiz1")
  } else if (type == "pbmc") {
    markerGenes <- c("Spi1", "Cd19", "Cd2", "Cd6", "Cd8a", "Cd4", "Csf1r", "S100a8", "Gata1", "Cd79a", "Cd3g", "Elane", "Mpo", "Itgam", "Retnlg",
                     "Ly6g", "S100a8", "Mpeg1", "Fn1", "Irf8", "Lyz2", "Cd74", "Rora", "Gzmk", "Cd3e", "Cd3d", "Meis1", "Kir2dl4", "Klrb1")
  }
  if (species == "human") {
    markerGenes <- toupper(markerGenes)
  } else if (species == "mouse") {
    markerGenes <- markerGenes
  }
}

########################################################################################################
#' @title regressCovBias
#' @description Calculate the residuals of a logistic regression of cell coverage vs. each irlba dimension
#'
#' @param reduction Name of dimensionality reduction to apply function over
#' @param method Calculate residuals using "lm" or "gam" model
#' @param obj Amethyst object containing irlba matrix to regress coverage bias from
#'
#' @return Returns the residuals
#' @export
#' @importFrom dplyr mutate select
#' @importFrom tibble column_to_rownames
#' @importFrom stats lm
#' @importFrom mgcv gam
#' @examples obj@reductions[["irlba_regressed"]] <- regressCovBias(obj)
regressCovBias <- function(
    obj,
    reduction = "irlba",
    method = "lm") {

  if (is.null(obj@reductions[[reduction]])) {
    stop("runIrlba step must be performed before regressing coverage bias.")
  }
  if (is.null(obj@metadata$cov)) {
    stop("Cell info must be available in the metadata slot before regressing coverage bias.")
  }

  regress <- merge(obj@metadata |>
                     dplyr::mutate(cov = log(cov)) |>
                     dplyr::select(cov), obj@reductions[[reduction]], by = 0) |>
    tibble::column_to_rownames(var = "Row.names")

  if (method == "lm") {
    result <- as.data.frame(apply(regress[c(2:ncol(regress))], 2, function(x) {
      lm_model <- stats::lm(x ~ regress$cov)
      stats::residuals(lm_model)
    }))
  } else if (method == "gam") {
    result <- as.data.frame(apply(regress[c(2:ncol(regress))], 2, function(x) {
      gam_model <- mgcv::gam(x ~ s(regress$cov))
      stats::residuals(gam_model)
    }))
  }
  return(result)
}

############################################################################################################################
#' @title runCluster
#' @description Determine cluster membership with Rphenograph https://github.com/JinmiaoChenLab/Rphenograph
#'
#' @param obj Amethyst object to perform clustering on
#' @param reduction Name of dimensionality reduction to calculate over
#' @param method Options are "louvain", which utilizes the Rphenograph https://github.com/JinmiaoChenLab/Rphenograph package; or "leiden", which utilizes the igraph package https://igraph.org/r/doc/cluster_leiden.html.
#' @param k integer; number of nearest neighbors
#' @param colname Character; name of column where results will be stored in metadata
#' @return Adds cluster membership to the metadata file of the amethyst object
#' @importFrom Rphenograph Rphenograph
#' @importFrom RANN nn2
#' @importFrom igraph membership
#' @importFrom tibble column_to_rownames
#' @importFrom leiden leiden
#' @importFrom rlang sym
#' @export
#' @examples
#' \dontrun{
#'   obj <- runCluster(obj = obj, k = 30, reduction = "irlba", method = "louvain", colname = "clusters_k30_irlba_louvain")
#' }
runCluster <- function(obj,
                       k = 50,
                       reduction = "irlba",
                       method = "louvain",
                       colname = "cluster_id") {

  if (method == "louvain") {
    clusters <- Rphenograph::Rphenograph(obj@reductions[[reduction]], k = {{k}})
    clusters <- do.call(rbind, Map(data.frame,
                                   cell_id = row.names(obj@reductions[[reduction]]),
                                   cluster_id = paste(igraph::membership(clusters[[2]]))))
  } else if (method == "leiden") { # https://github.com/TomKellyGenetics/leiden

    cell_names <- rownames(obj@reductions[[reduction]])
    snn <- RANN::nn2(obj@reductions[[reduction]], k = {{k}})$nn.idx
    adjacency_matrix <- matrix(0L, length(cell_names), length(cell_names))
    rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- cell_names

    for(ii in 1:length(cell_names)) {
      adjacency_matrix[ii,cell_names[snn[ii,]]] <- 1L
    }
    clusters <- data.frame(
      row.names = cell_names,
      cell_id = cell_names,
      cluster_id = leiden::leiden(adjacency_matrix)
    )
    clusters$cluster_id <- as.character(clusters$cluster_id)

  } else {
    stop("Please choose method 'louvain' or 'leiden'. Any alternative method may be implemented manually and added to the object metadata.")
  }

  if (colname != "cluster_id") {
    clusters <- clusters |> dplyr::rename(!!sym(colname) := "cluster_id")
  }

  if (is.null(obj@metadata[[colname]])) {
    if (is.null(obj@metadata)) {
      obj@metadata <- clusters |> dplyr::select(-cell_id)
    }
    if (!is.null(obj@metadata)) {
      obj@metadata <- merge(obj@metadata, clusters, by = 0) |> tibble::column_to_rownames(var = "Row.names")
      obj@metadata <- obj@metadata[, !grepl("cell_id", colnames(obj@metadata))]
    }
  } else if (!is.null(obj@metadata[[colname]])) {
    obj@metadata <- merge(obj@metadata |> dplyr::select(-!!sym(colname)),
                          clusters, by = 0) |>
      tibble::column_to_rownames(var = "Row.names")
    obj@metadata <- obj@metadata[, !grepl("cell_id", colnames(obj@metadata))]
  }
  output <- obj
}
