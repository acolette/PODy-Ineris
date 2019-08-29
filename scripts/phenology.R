# load need package
library(ncdf4)
library(lubridate)

# some custom functions to read the input data
source('scripts/functions.R')

# define input/output filename & select the specie
specie <- "Wheat" #c("Wheat", "Tomato", "Potato", "Grass", "Beech", "TemperateOak", "Spruce") 
ifile  <- "data/EDT/pody_input_EDT_CHIF_2010.nc" # met soil pol data
dfile  <- "data/INPUT/DEPO_PARS"                 # depot parameter
sfile  <- "data/INPUT/SPECIES_PAR.csv"           # specie parameter
ofile  <- paste0("output/phenology_EDT_2010_", specie, ".dat")

# read input data ---------------------------------------------------------
# read phenology & deposition parameters
p_specie <- read_param_specie(sfile, specie)
p_depo   <- read_param_deposition(dfile, specie)

# read data (& day by day temperature in the loop)
fic <- nc_open(ifile)
# dimensions (x,y,t)
lon   <- ncvar_get(fic, "lon")
lat   <- ncvar_get(fic, "lat")
times <- ncvar_get(fic, "Times") 
times <- ymd_hms(times, tz = "UTC") # chr => POSIXct
nday  <- length(unique(date(times)))
nlon  <- length(lon[,1])
nlat  <- length(lon[1,])
# soil property
alti <- ncvar_get(fic, "z")  # alitude (m)
PWP  <- ncvar_get(fic, "WP") # Permament Wilting Point ([0-1])
FC   <- ncvar_get(fic, "FC") # Field Capacity ([0-1])

# phenology computation ---------------------------------------------------
print(paste0("Compute phenology for ", specie))

# Variable declaration & some computation required once
fphen     <- array(0, dim = c(nlon, nlat, nday))
tem2_day  <- array(0, dim = c(nlon, nlat, nday))
tmp_fphen <- array(0, dim = c(nlon, nlat))
dlai      <- numeric(nday) # unique value for all the domain, could be gridded 
dsai      <- numeric(nday) # unique value for all the domain, could be gridded 

if (specie %in% c("Wheat", "Potato", "Tomato")) {
  # init origin of degdays:  midanthesis in degdays(Wheat), 
  #  TuberInialization (Potato) & TomatoPlantation  (tomato) fixed days
  midanthesis2d <- array(1075, dim = c(nlon, nlat))
  # fixed day at the 1st of june (cf Mapping Maual) 
  Potatoday <- 151 # unique value for all the domain, could be gridded 
  Tomatoday <- 151 # unique value for all the domain, could be gridded 
  
} else if (specie  %in% c("Beech", "TemperateOak")) {
  # latitude model for TREES for EGS and SGS see mapping manual PODySPEC
  Astart_FD <- array(0, dim = c(nlon, nlat))
  Aend_FD   <- array(0, dim = c(nlon, nlat))
  # Astart_FD and Aend_FD are given from the EMEP latitude model for forest trees 
  # and corrected with the altitude parameter according to the MM
  # Effect of altitude on phenology is incorporated by assuming a later SGS 
  # & earlier EGS by 10 days for every 1000 m
  Astart_FD <- 105 - 2 * (50 - lat) + 10 * round(alti / 1000, 0)
  Aend_FD   <- 297 - 2 * (lat - 50) - 10 * round(alti / 1000, 0)
} 

# loop over day
for (d in seq_len(nday)) {
  # read temperature (hourly, degC) avg to daily & compute degdays
  ind_start   <- c(1, 1, (d - 1) * 24 + 1) # ind to read only 24h 
  tem2_hour   <- ncvar_get(fic, "tem2", start = ind_start, count = c(nlon,nlat,24))
  tem2_oneday <- rowMeans(tem2_hour, dims = 2) # dayavg
  # degdays with base temp = 0
  tem2_oneday[which(tem2_oneday < 0)] <- 0 
  tem2_day[,,d] <- tem2_oneday
  # compute degdays sum 
  cumtem <- apply(tem2_day[,,1:d], 1:2, sum) 
  
  # computation of lai and sai according to chimere methodolology
  if (d >= p_depo$depsgs) {
    dlai[d] <- p_depo$deplai2
    dsai[d] <- dlai[d] + 1.5
    
    if (d <= (p_depo$depsgs + p_depo$depsgl)) {
      dlai[d] <- p_depo$deplai1 + (d - p_depo$depsgs) * 
        (p_depo$deplai2 - p_depo$deplai1) / p_depo$depsgl
      dsai[d] <- 5 * dlai[d] / 3.5
    }
    
    if (specie %in% c("Spruce","Beech", "TemperateOak")) { 
      dsai[d] <- dlai[d] + 1
    } #trees
 
    if (specie == "Grass") {
      dsai[d] <- dlai[d]
    }
  } # lai & sai

  # CROP --------------------------------------------------------------------
  # CROP inside -> phenology function for Wheat, Potato and Tomato
  if (specie %in% c("Wheat", "Potato", "Tomato")) {
    # init specie cumtem variable 
    PotatoCumtem  <- array(0, dim = c(nlon, nlat))
    TomatoDegDay  <- array(0, dim = c(nlon, nlat))
    
    # Computation of degredays for phenology application and accumulation period. 
    # Origin is midanthesis or TuberInialization or TomatoPlantation
    if (specie == "Wheat") {
      cumtem <- cumtem - midanthesis2d #degdays from mid anthesys
    }
    if (specie == "Potato") {
      PotatoCumtem <- apply(tem2_day[,, 1:Potatoday], 1:2, sum)
    }
    if (specie == "Tomato") {
      if (d < Tomatoday) {
        # For day below Tomatoday, no phenology -> Astart_ETS -1 ==> no phen loop
        TomatoDegDay <- array(p_specie$Astart_ETS - 1, dim = c(nlon, nlat)) 
      } else {
        TomatoDegDay <- tem2_day[,, Tomatoday:d]
        TomatoDegDay <- TomatoDegDay - 10 # Base temperature 10°C
        TomatoDegDay[which(TomatoDegDay < 0)] <- 0
        TomatoDegDay <- apply(TomatoDegDay, c(1, 2), sum) # degdays(base temperature = 10°C) from Tomatoday
      }
    } # tomato
    
    # phenology computation (crops) -------------------------------------------
    if (specie == "Wheat") {
      tmp_fphen <- fphen[,,d]
      # condition n°1
      this_whichWheatN1 <- which((cumtem > (p_specie$fphen_2_ets + p_specie$fphen_1_ets)) & 
                                   (cumtem < (p_specie$fphen_2_ets + p_specie$fphen_3_ets))) 
      
      if (length(this_whichWheatN1) > 0) {
        tmp_fphen[this_whichWheatN1] <- 1
      }
      # condition n°2
      this_whichWheatN2 <- which((cumtem > (p_specie$fphen_2_ets + p_specie$fphen_3_ets)) & 
                                   (cumtem < (p_specie$fphen_2_ets + p_specie$fphen_4_ets)))
      if (length(this_whichWheatN2) > 0) {
        tmp_fphen[this_whichWheatN2] <- 1 - (p_specie$fphen_a / (p_specie$fphen_4_ets - p_specie$fphen_3_ets)) * 
          (cumtem[this_whichWheatN2] - p_specie$fphen_3_ets)
      }
      # condition n°3
      this_whichWheatN3 <- which((cumtem > (p_specie$fphen_2_ets + p_specie$fphen_4_ets)) & 
                                   (cumtem < p_specie$fphen_5_ets))
      if (length(this_whichWheatN3) > 0) {
        tmp_fphen[this_whichWheatN3] <- p_specie$fphen_e - (p_specie$fphen_e / 
          (p_specie$fphen_5_ets - p_specie$fphen_4_ets) * (cumtem[this_whichWheatN3] - p_specie$fphen_4_ets))
      }
      
      if (length(tmp_fphen) != 0) {
        fphen[,,d] <- tmp_fphen
      }
    } # Wheat
  
    if (specie == "Potato") {
      tmp_fphen <- fphen[,,d]
      # fphen from Background Document A, Nov 2017 p27, Fig 4.2 (Pleijel et al., 2007)
      # condition n°1
      this_wichPotatoN1 <- which(((cumtem - PotatoCumtem) >= p_specie$fphen_1_ets) & ((cumtem - PotatoCumtem) < 0))
      if (length(this_wichPotatoN1) > 0) {
        tmp_fphen[this_wichPotatoN1] <- 1 - (((1 - p_specie$fphen_a) / p_specie$fphen_1_ets) * 
                                               (cumtem[this_wichPotatoN1] - PotatoCumtem[this_wichPotatoN1]))
      } 

      # condition n°2
      this_wichPotatoN2 <- which(((cumtem - PotatoCumtem) >= 0) & ((cumtem - PotatoCumtem) <= p_specie$fphen_2_ets))
      if (length(this_wichPotatoN2) > 0) {
        tmp_fphen[this_wichPotatoN2] <- 1 - ((1 - p_specie$fphen_e) / p_specie$fphen_2_ets) * 
          (cumtem[this_wichPotatoN2] - PotatoCumtem[this_wichPotatoN2])
      }
      
      if (length(tmp_fphen) != 0) {
        fphen[,,d] <- tmp_fphen
      }
    } # Potato
    
    if (specie == "Tomato") {
      tmp_fphen <- fphen[,,d]
      # fphen from Background Document A, Nov 2017 p29, Fig 4.7 (González-Fernández et al., 2014)
      # condition n°1
      this_wichTomatoN1 <- which((TomatoDegDay >= p_specie$Astart_ETS) & (TomatoDegDay <= p_specie$Aend_ETS))
      if (length(this_wichTomatoN1) > 0) {
        tmp_fphen[this_wichTomatoN1] <- (TomatoDegDay[this_wichTomatoN1] - p_specie$fphen_2_ets) / 
          (p_specie$Astart_ETS - p_specie$fphen_2_ets)
      }
      
      if (length(tmp_fphen) != 0) {
        fphen[,,d] <- tmp_fphen
      }
    } # tomato 
  } # if crops
  
  # phenology (tree = fixed day) --------------------------------------------
  if (specie  %in% c("Beech", "TemperateOak")) {

    tmp_fphen <- fphen[,,d]
    this_whichBeechN1 <- which(Astart_FD >= d)
    if (length(this_whichBeechN1) > 0) {
      tmp_fphen[this_whichBeechN1] <- p_specie$fphen_a
    }
    
    this_whichBeechN2 <- which((Astart_FD < d) & (Astart_FD >= (d - p_specie$fphen_1_FD)))
    if (length(this_whichBeechN2) > 0) {
      tmp_fphen[this_whichBeechN2] <- ((1 - p_specie$fphen_a) * ((d - Astart_FD[this_whichBeechN2]) / 
                                                                   p_specie$fphen_1_FD) + p_specie$fphen_a)
    }
    
    if (d <= p_specie$LIMstart_FD) {
      this_whichBeechN3 <- which((Astart_FD < (d - p_specie$fphen_1_FD)))
    } else {
      this_whichBeechN3 <- NULL
    }
    
    if (length(this_whichBeechN3) > 0) {
      tmp_fphen[this_whichBeechN3] <- p_specie$fphen_b
    }
    
    if ((d > p_specie$LIMstart_FD) & (d < (p_specie$LIMstart_FD + p_specie$fphen_2_FD))) {
      tmp_fphen[] <- (1 - p_specie$fphen_c) * (((p_specie$fphen_2_FD + p_specie$LIMstart_FD) - 
        (p_specie$LIMstart_FD + (d - p_specie$LIMstart_FD))) / p_specie$fphen_2_FD) + p_specie$fphen_c
    }
    
    if ((d >= (p_specie$LIMstart_FD + p_specie$fphen_2_FD)) & (d <= (p_specie$LIMend_FD - p_specie$fphen_3_FD))) {
      tmp_fphen[] <- p_specie$fphen_c
    }
    
    if ((d > p_specie$LIMend_FD - p_specie$fphen_3_FD) & (d <= p_specie$LIMend_FD)) {
      tmp_fphen[] <- (1 - p_specie$fphen_c) * ((d - (p_specie$LIMend_FD - 
        p_specie$fphen_3_FD)) / p_specie$fphen_3_FD) + p_specie$fphen_c
    }
    
    if (d >= p_specie$LIMend_FD) {
      this_whichBeechN4 <- which(Aend_FD >= d + p_specie$fphen_4_FD)
    }
    
    if (length(this_whichBeechN4) > 0) {
      tmp_fphen[this_whichBeechN4] <- p_specie$fphen_d
    }

    this_whichBeechN5 <- which((Aend_FD > d) & (Aend_FD < (d - p_specie$fphen_4_FD)))
    if (length(this_whichBeechN5) > 0) {
      tmp_fphen[this_whichBeechN5] <- (1 - p_specie$fphen_e) * 
        ((Aend_FD[this_whichBeechN5] - d) / p_specie$fphen_4_FD) + p_specie$fphen_e
    }
    
    this_whichBeechN6 <- which(Aend_FD <= d)
    if (length(this_whichBeechN6) > 0) {
      tmp_fphen[this_whichBeechN6] <- p_specie$fphen_e
    }

    if (length(tmp_fphen) != 0) {
      fphen[,,d] <- tmp_fphen
    }
  } # if trees 
  
  # Norway spruce in the continental region, the growth period is determined by 
  # air temperature defined according to the ftemp function.
  # The growing season is assumed to occur when air temperatures are between the 
  # Tmin and Tmax thresholds of the ftemp relationships. 
  # During such periods there is no limitation on conductance associated with 
  # leaf development stage (i.e. fphen = 1)
  if (specie %in% c("Spruce", "Grass")) {
    fphen[,,d] <- 1
  }
} # day loop

# save phenology & some other inputs/informations
# common variables (all specie)
common_var <- c("lon", "lat", "nlon", "nlat", "nday", "times", "PWP", "FC", 
                "tem2_day",  "dsai", "dlai", "p_depo", "p_specie", "fphen")
# species related variables
specie_var <- list(
  "Wheat"        = "midanthesis2d",
  "Tomato"       = "Tomatoday",
  "Potato"       = "Potatoday",
  "Beech"        = c("Astart_FD", "Aend_FD"),
  "TemperateOak" = c("Astart_FD", "Aend_FD")
)
# save all variables
save(list = c(common_var, specie_var[[specie]]), file = ofile)

# close nc file
nc_close(fic) 