############################################################################################################################
#' @title createObject
#' @description Create object of class amethyst
#'
#' @param h5paths Path to the hdf5 file containing base-level read information organized by methylation type and barcode
#' @param index Corresponding gene coordinates in the hdf5 file
#' @param metadata Optional cell metadata. If included, make sure row names are cell IDs
#' @param ref Genome annotation file with chromosome, start, and end position information for genes of interest. See the "makeRef" function.
#' @return Returns a single object of class amethyst with h5path, index, metadata, ref, reduction, and irlba slots.
#' @importFrom methods new
#' @export
#' @examples obj <- createObject(h5paths = "~/Downloads/test.h5")
createObject <- function(h5paths = NULL,
                         index = NULL,
                         ref = NULL,
                         metadata = NULL) {
  methods::new(Class = "amethyst")
}

methods::setClass("amethyst", slots = c(
  h5paths = "ANY",
  genomeMatrices = "ANY",
  irlba = "ANY",
  index = "ANY",
  metadata = "ANY",
  ref = "ANY"
))

############################################################################################################################
#' @title extractAttributes
#' @description Extract information from the attributes column of a gtf file.
#' This helper function was taken from https://www.biostars.org/p/272889/
#'
#' @param gtf_attributes A string containing the attributes column from a GTF file
#' @param att_of_interest Name of the attribute containing the desired information
#' @return Returns the value associated with the specified attribute key from the GTF attributes string
#' @export
#' @examples gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))
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
#' @description Generate an annotation file from a gtf. Paths to mm10 and hg38 are hard-coded.
#' Please provide the gtf if analyzing other species.
#'
#' @param ref Reference genome, i.e. "hg38" or "mm10"
#' @param gtf If not using hg38 or mm10, please provide the gtf.gz file path
#' @param attributes Information to extract from the gtf file. Must be a column name
#' @return Returns an annotated reference of gene locations
#' @importFrom dplyr mutate
#' @importFrom rtracklayer readGFF
#' @export
#' @examples ref <- makeRef(ref = "hg38")
makeRef <- function(ref,
                    gtf = NULL,
                    attributes = c("gene_name", "exon_number")) {
  if (ref == "hg38") {
    options(timeout = 1000)
    gtf <- rtracklayer::readGFF("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz")
  } else if (ref == "mm10") {
    options(timeout = 1000)
    gtf <- rtracklayer::readGFF("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz")
  } else if (gtf) {
    gtf <- rtracklayer::readGFF(gtf)
  } else {
    print("Ref was different from default (hg38) but not provided. Please provide unzipped gtf reference annotation file.")
  }
  for (i in attributes) {
    gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))
  }
  gtf <- gtf |>
    dplyr::mutate(location = paste0(seqid, "_", start, "_", end))
}

############################################################################################################################
#' @title fetchMarkers
#' @description Fetch a pre-selected list of key marker genes
#' Key marker genes for brain methylation selected from http://neomorph.salk.edu/omb/ct_table
#'
#' @param ref If using "hg38" or "mm10"
#' @param type Tissue type of interest. For now, "brain" and "pbcm" are options.
#' @return Returns a character vector of key marker genes.
#' @export
#' @examples markerGenes <- fetchMarkers(ref = "hg38", type = "brain")
fetchMarkers <- function(ref,
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

  }
  if (ref == "hg38") {
    markerGenes <- toupper(markerGenes)
  } else if (ref == "mm10") {
    markerGenes <- markerGenes
  }
}

############################################################################################################################
#' @title dimEstimate
#' @description Estimate the nv value needed for singular value decomposition with irlba
#'
#' @param obj Amethyst object containing the matrix to calculate, which should be in the genomeMatrices slot
#' @param genomeMatrices Name of the matrix in the genomeMatrices slot to calculate nv
#' @param threshold Amount of variance that must be explained by the nv value, with 1 being 100% of variance.
#' @param dims Number of singular values to test for each matrix
#'
#' @return Integer indicating the number of principal components required to meet the variance threshold for each matrix
#' @export
#'
#' @examples dimEstimate(obj = combined, genomeMatrices = c("cg_100k_score", "ch_100k_pct"), dims = c(50, 50), threshold = 0.98)
#'
dimEstimate <- function(
    obj,
    genomeMatrices,
    dims,
    threshold = 0.98) {

  if (length(genomeMatrices) != length(dims)) {
    stop("Number of input matrices must equal the length of the dimension list")
  }

  svd_output <- list()
  for (i in 1:length(genomeMatrices)) {
    svd_output[[i]] <- obj@genomeMatrices[[genomeMatrices[i]]]
    windows <- rownames(svd_output[[i]])
    cells <- colnames(svd_output[[i]])
    svd_output[[i]][is.na(svd_output[[i]])] <- 0
    svd_output[[i]] <- irlba::irlba(as.matrix(svd_output[[i]]), dims[[i]])
    svd_output[[i]] <- as.data.frame(svd_output[[i]]$d)
    colnames(svd_output[[i]]) <- paste(genomeMatrices[[i]])
    rownames(svd_output[[i]]) <- paste("DIM", 1:dims[[i]])
  }
  svd_output <- do.call(cbind, svd_output)

  dims_to_use <- apply(svd_output, 2, function(x) {
    cumulative_variance_explained <- cumsum((x^2) / sum(x^2))
    if (any(cumulative_variance_explained >= threshold)) {
      min(which(cumulative_variance_explained >= threshold))
    } else {
      stop(paste("No nv value explains ", (threshold*100), "% of variance."))
    }
  })

  dims_to_use

}

############################################################################################################################
#' @title runIrlba
#' @description Perform dimensionality reduction based on methylation levels over a matrix stored in the @genomeMatrices slot
#'
#' @param obj Object for which to run irlba
#' @param genomeMatrices list of matrices in the genomeMatrices slot to use for irlba
#' @param dims list of how many dimensions to output for each matrix
#' @return Returns a matrix of appended irlba dimensions as columns and cells as rows
#' @importFrom irlba irlba
#' @export
#' @examples obj <- runIrlba(obj, genomeMatrices = c("ch_2M_pct", "ch_2M_score", "cg_2M_score"), dims = c(10, 10, 10))
runIrlba <- function(
    obj,
    genomeMatrices,
    dims) {

  if (length(genomeMatrices) != length(dims)) {
    stop("Number of input matrices must equal the length of the dimension list")

  } else {
    if (!is.null(obj@irlba)) {
      obj@irlba <- NULL
    }

    matrix <- list()

    for (i in 1:length(genomeMatrices)) {
      matrix[[i]] <- obj@genomeMatrices[[genomeMatrices[i]]]
      windows <- rownames(matrix[[i]])
      cells <- colnames(matrix[[i]])
      matrix[[i]][is.na(matrix[[i]])] <- 0
      matrix[[i]] <- irlba::irlba(as.matrix(matrix[[i]]), dims[[i]])
      matrix[[i]] <- as.data.frame(matrix[[i]]$v)
      colnames(matrix[[i]]) <- paste("DIM", 1:dims[[i]], "_", genomeMatrices[i], sep = "")
      rownames(matrix[[i]]) <- cells
    }

    obj@irlba <- do.call(cbind, matrix)
    return(obj)
  }
}

########################################################################################################
#' @title regressCovBias
#' @description Calculate the residuals of a logistic regression of cell coverage vs. each irlba dimension
#'
#' @param obj Amethyst object containing irlba matrix to regress coverage bias from
#' @return Returns the residuals in the irlba slot
#' @export
#' @importFrom dplyr mutate select
#' @importFrom tibble column_to_rownames
#' @importFrom stats lm
#' @examples obj <- regressCovBias(obj)
regressCovBias <- function(
    obj) {

  if (is.null(obj@irlba)) {
    stop("runIrlba step must be performed before regressing coverage bias.")
  }
  if (is.null(obj@metadata$cov)) {
    stop("Cell info must be available in the metadata slot before regressing coverage bias.")
  }

  regress <- merge(obj@metadata |> dplyr::mutate(cov = log(cov)) |> dplyr::select(cov), obj@irlba, by = 0) |> tibble::column_to_rownames(var = "Row.names")
  obj@irlba <- as.data.frame(apply(regress[c(2:ncol(regress))], 2, function(x) {
    lm_model <- stats::lm(x ~ regress$cov)
    stats::residuals(lm_model)
  }))

  return(obj)
}

############################################################################################################################
#' @title runCluster
#' @description Determine cluster membership with Rphenograph https://github.com/JinmiaoChenLab/Rphenograph
#'
#' @param obj Amethyst object to perform clustering on
#' @param k_phenograph integer; number of nearest neighbors
#' @return Adds cluster membership to the metadata file of the amethyst object
#' @importFrom Rphenograph Rphenograph
#' @importFrom igraph membership
#' @importFrom tibble column_to_rownames
#' @export
#' @examples obj <- runCluster(obj = obj, k_phenograph = 30)
runCluster <- function(obj,
                       k_phenograph = 50) {

  clusters <- Rphenograph::Rphenograph(obj@irlba, k = k_phenograph)
  clusters <- do.call(rbind, Map(data.frame, cell_id = row.names(obj@irlba), cluster_id = paste(igraph::membership(clusters[[2]]))))

  if (is.null(obj@metadata$cluster_id)) {
    if (is.null(obj@metadata)) {
      obj@metadata <- clusters |> dplyr::select(-cell_id)
    }
    if (!is.null(obj@metadata)) {
      obj@metadata <- merge(obj@metadata |> dplyr::select(-cluster_id), clusters, by = 0) |> tibble::column_to_rownames(var = "Row.names")
      obj@metadata <- obj@metadata[, !grepl("cell_id", colnames(obj@metadata))]
    }
  } else if (!is.null(obj@metadata$cluster_id)) {
    obj@metadata <- merge(obj@metadata |> dplyr::select(-cluster_id), clusters, by = 0) |> tibble::column_to_rownames(var = "Row.names")
    obj@metadata <- obj@metadata[, !grepl("cell_id", colnames(obj@metadata))]
  }
  output <- obj
}


############################################################################################################################
#' @title runUmap
#' @description Perform dimension reduction with Uniform Manifold APproximation and Projection for Dimension Reduction
#'
#' @param obj Amethyst object for which to determine umap coordinates
#' @param neighbors Number of closest points to factor into projection calculations. A higher number will capture more global structure
#' @param dist Distance between point pairs
#' @param method Distance metric to utilize. Default is euclidean
#' @return Adds umap_x and umap_y coordinates to the metadata file of the object
#' @importFrom umap umap
#' @export
#' @examples obj <- runUmap(obj, neighbors = 30, dist = 0.1, method = "euclidean")
runUmap <- function(obj,
                    neighbors = 30,
                    dist = 0.1,
                    method = "euclidean") {

  umap_dims <- as.data.frame(umap::umap(as.data.frame(obj@irlba), method = "naive", dims = 2, n_components = 2, neigh = neighbors, mdist = dist, metric = method)$layout)

  # Add to metadata
  if (is.null(obj@metadata$umap_x)) {
    metadata <- merge(umap_dims, obj@metadata, by = 0) %>%
      tibble::column_to_rownames(var = "Row.names") %>%
      dplyr::rename("umap_x" = "V1", "umap_y" = "V2")
  } else {
    metadata <- merge(umap_dims, obj@metadata %>% dplyr::select(-c(umap_x, umap_y)), by = 0) %>%
      tibble::column_to_rownames(var = "Row.names") %>%
      dplyr::rename("umap_x" = "V1", "umap_y" = "V2")
  }
  obj@metadata <- metadata
  output <- obj
}
