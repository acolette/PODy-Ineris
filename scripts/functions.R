# read data ---------------------------------------------------------------
read_param_specie <- function(ifile, specie) {
  # read parameter file
  # define deposition parameters (include growing season for each vegetation type)
  # obtained from Chimere & adapted with WGE Mapping Manual
  specie_par <- read.table(ifile, sep = ";", header = TRUE, stringsAsFactors = FALSE)
  specie_par <- specie_par[,-c(1,2)]
  specie_par <- apply(specie_par, 2, as.numeric)

  # test presence absence of the specie
  specie_col <- colnames(specie_par)
  if (!specie %in% specie_col) {
    stop("specie not define in input parameters file")
  } 
  
  # return all args as list
  list(
    gmax          = specie_par[ 1, specie],
    fmin          = specie_par[ 2, specie],
    light_a       = specie_par[ 3, specie],
    Tmin          = specie_par[ 4, specie],
    Topt          = specie_par[ 5, specie],
    Tmax          = specie_par[ 6, specie],
    VPDmax        = specie_par[ 7, specie],
    VPDmin        = specie_par[ 8, specie],
    VPDcrit       = specie_par[ 9, specie], 
    # PAWt used for all crops (not used for trees)
    # PAWt is rescaled between [0-1] 
    PAWt          = specie_par[10, "Wheat"] / 100, 
    fO3cst        = specie_par[15, specie],
    fO3AOT        = specie_par[16, specie],
    fO3exp        = specie_par[17, specie],
    Astart_ETS    = specie_par[18, specie],
    Aend_ETS      = specie_par[19, specie],
    leafdimension = specie_par[20, specie],
    zcanopy       = specie_par[21, specie],
    fphen_a       = specie_par[22, specie],
    fphen_b       = specie_par[23, specie],
    fphen_c       = specie_par[24, specie],
    fphen_d       = specie_par[25, specie],
    fphen_e       = specie_par[26, specie],
    fphen_1_ets   = specie_par[27, specie],
    fphen_2_ets   = specie_par[28, specie],
    fphen_3_ets   = specie_par[29, specie],
    fphen_4_ets   = specie_par[30, specie],
    fphen_5_ets   = specie_par[31, specie],
    fphen_1_FD    = specie_par[32, specie],
    fphen_2_FD    = specie_par[33, specie],
    fphen_3_FD    = specie_par[34, specie],
    fphen_4_FD    = specie_par[35, specie],
    LIMstart_FD   = specie_par[36, specie],
    LIMend_FD     = specie_par[37, specie]
  )
}
read_param_deposition <- function(ifile, specie) {
  # define corresp chimere column and specie list
  # use for computing SAI & LAI according to Chimere methodology
  chimcol <- c("Wheat" = 17, "Potato" = 18, "Tomato" = 20, "Spruce" = 23, 
               "Beech" = 21, "TemperateOak" = 21, "Grass" = 19)
  
  # test presence absence of the specie
  if (!specie %in% names(chimcol)) {
    stop("specie not define in chimere column")
  } 
  
  # read param file & return all args as list
  depo_pars <- read.table(ifile, header = FALSE, skip = 1)
  list(
    depsgs  = depo_pars[ 9, chimcol[specie]],
    depegs  = depo_pars[10, chimcol[specie]],
    depsgl  = depo_pars[11, chimcol[specie]],
    depegl  = depo_pars[12, chimcol[specie]],
    deplai1 = depo_pars[13, chimcol[specie]],
    deplai2 = depo_pars[14, chimcol[specie]],
    depphe0 = depo_pars[15, chimcol[specie]],
    depphe1 = depo_pars[16, chimcol[specie]],
    depphe2 = depo_pars[17, chimcol[specie]]
  )
}




