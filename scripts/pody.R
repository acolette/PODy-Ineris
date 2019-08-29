library(ncdf4)
library(lubridate)

# user define args & input/output files
specie <- "Spruce" # "TemperateOak" # "Beech" #"Grass" #"Potato" # "Tomato" # Wheat" # which species
specie_list  <- c("Wheat", "Tomato", "Potato", "Grass", "Beech", "TemperateOak", "Spruce") 
O3layerhght  <- 10 # O3 layer height (downscalling if layer > canopy)
heavy_output <- TRUE # netcdf output if TRUE 

ifile  <- "data/EDT/pody_input_EDT_CHIF_2010.nc"                # netcdf
pfile  <- paste0("output/phenology_EDT_2010_", specie, ".dat")  # phenology 
ofile  <- paste0("output/pody_input_EDT_2010_", specie, ".dat") # ofile 

# init variables ----------------------------------------------------------
# constant
avogadro <- 6.022045e23

# appropriate Y threeshold for species 
if (specie %in% c("Wheat", "Potato", "Tomato")) {
  y_pod <- 6
} else {
  # Trees and Grass
  y_pod <- 1
}

# read data --------------------------------------------------------------
# phenology (depends of year and specie)
load(pfile)

# init array
# POD and AOT
pody   <- array(0, dim = c(nlon, nlat)) # used to store final result
pody0  <- array(0, dim = c(nlon, nlat)) # used for fO3 computation
AOT0   <- array(0, dim = c(nlon, nlat)) # used for fO3 computation

# Gmax limitation functions
ftemp  <- array(0, dim = c(nlon, nlat)) # ftemp
fO3    <- array(1, dim = c(nlon, nlat)) # fO3
fsmi   <- array(1, dim = c(nlon, nlat)) # soil moisture index   
# Resistances parameter
Rext   <- array(2500, dim = c(nlon, nlat))
Rsoil  <- array(200, dim = c(nlon, nlat))
O3Gext <- 1 / Rext # mapping manual methodology (conductance ext)

# netcdf file 
fic <- nc_open(ifile)

# compute pody --------------------------------------------------------------
print(paste0("Compute POD", y_pod, specie))
dh <- 0  # init hour 

# netcdf declaration ------------------------------------------------------
if (heavy_output) {
  # define dimensions
  londim  <- ncdim_def("lon", "degrees_east", lon[,1])
  latdim  <- ncdim_def("lat", "degrees_north", lat[1,])
  timedim <- ncdim_def("time", paste0("hour nb from ", unique(year(times))), seq_len(length(times)))
  chardim <- ncdim_def("DateStrLen", "", 1:19, create_dimvar = FALSE) # to store Times axis
  
  # define variables
  fillvalue <- NA
  # 1D variables
  def_tim <- ncvar_def(name = "Times", units = "", dim = list(chardim, timedim),prec = "char")
  # 2D variables
  def_lon <- ncvar_def(name = "lon_2d", units = "deg", dim = list(londim,latdim),
               missval = fillvalue, longname = "longitude 2d", prec = "single")
  def_lat <- ncvar_def(name = "lat_2d", units = "deg", dim = list(londim,latdim),
               missval = fillvalue, longname = "latitude 2d", prec = "single")
  # 3D variables
  def_pody   <- ncvar_def(name = "PODy", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = paste0("POD", y_pod, specie), prec = "single")
  def_O3     <- ncvar_def(name = "O3", units = "ppb", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Ozone canopy", prec = "single")
  def_fO3    <- ncvar_def(name = "fO3", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Ozone limitation function", prec = "single")
  def_fsmi   <- ncvar_def(name = "fsmi", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Soil moisture limitation function", prec = "single")
  def_ftemp  <- ncvar_def(name = "ftemp", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Temperature limitation function", prec = "single")
  def_flig   <- ncvar_def(name = "flig", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Light limitation function", prec = "single")
  def_fphen  <- ncvar_def(name = "fphen", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Phenology limitation function", prec = "single")
  def_fvpd   <- ncvar_def(name = "fvpd", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Vapour pressure limitation function", prec = "single")
  def_accper <- ncvar_def(name = "accper", units = "", dim = list(londim,latdim,timedim),
                  missval = fillvalue, longname = "Accumulation period", prec = "single")
  
  
  # create netCDF file and put arrays
  ncfname <- paste0("output/pody_input_EDT_2010_", specie, ".nc")
  ncout   <- nc_create(ncfname, list(def_tim, def_pody,def_O3, def_fO3, def_fsmi, 
                                     def_ftemp, def_flig,def_fphen, def_fvpd, 
                                     def_accper), force_v4 = TRUE)
}

for (d in seq_len(nday)) {
  print(d)
  
  # init variables
  accperiod <- array(0, dim = c(nlon, nlat)) # used to save period
  
  # recompute cumtem for day d
  cumtem <- apply(tem2_day[,, 1:d], 1:2, sum)

  # Now Time origin is: midanthesis (wheat), init of the tubercule (potato) & 
  # transplantation into the field (tomato)
  if (specie == "Wheat") {
    cumtem <- cumtem - midanthesis2d # degdays from mid anthesys
  }
  if (specie == "Potato") {
    PotatoCumtem <- apply(tem2_day[,, 1:Potatoday], 1:2, sum)
  }
  if (specie == "Tomato") {
    if (d < Tomatoday) {
      # For day below Tomatoday, TomatoDegDay <- Astart_ETS -1 => never  in the accumulation period
      TomatoDegDay <- array(p_specie$Astart_ETS - 1, dim = c(nlon, nlat))
    } else {
      # accumulation period for tomato is 250ºC days to 1500ºC days after 
      # transplantation in the field over a base temperature of 10ºC (González-Fernández et al., 2014)
      TomatoDegDay <- tem2_day[,, Tomatoday:d]
      TomatoDegDay <- TomatoDegDay - 10 # Base temperature 10°C
      TomatoDegDay[which(TomatoDegDay < 0)] <- 0
      TomatoDegDay <- apply(TomatoDegDay, c(1,2), sum) # degdays(base temp = 10°C) from Tomatoday
    }
  }
  
  # loop over hours
  for (h in 0:23) {
    # init & count hour
    this_which <- NULL # reinit this_which
    dh <- dh + 1
    
    # reinitilize the limitation functions 
    fsmi  <- array(1, dim = c(nlon, nlat))
    ftemp <- array(1, dim = c(nlon, nlat))
    flig  <- array(1, dim = c(nlon, nlat))
    fvpd  <- array(1, dim = c(nlon, nlat))
    
    # get input variables -------------------------------------------------
    nc_ind  <- list(start = c(1, 1, dh), count = c(nlon, nlat, 1))
    O3hppb  <- ncvar_get(fic,"O3",    start = nc_ind[["start"]], count = nc_ind[["count"]]) # ppb needed
    swvl    <- ncvar_get(fic, "swvl2",start = nc_ind[["start"]], count = nc_ind[["count"]]) # [0-1]
    ustaloc <- ncvar_get(fic, "usta", start = nc_ind[["start"]], count = nc_ind[["count"]]) # m/s
    srehloc <- ncvar_get(fic, "sreh", start = nc_ind[["start"]], count = nc_ind[["count"]]) # [0,1]
    tem2loc <- ncvar_get(fic, "tem2", start = nc_ind[["start"]], count = nc_ind[["count"]]) # in degC
    airmloc <- ncvar_get(fic, "airm", start = nc_ind[["start"]], count = nc_ind[["count"]]) # in molec/cm3
    obukloc <- ncvar_get(fic, "obuk", start = nc_ind[["start"]], count = nc_ind[["count"]]) # Obukhov length (m)
    swrd    <- ncvar_get(fic, "swrd", start = nc_ind[["start"]], count = nc_ind[["count"]]) # swrd W.m-2
    ppfd    <- swrd * 4.5 * 0.5 # (Photosynthetic Photon Flux Density, umol/m2/s)
    # Here we assume (1) 4.5 umol/m2/s per W/m2 (changed in MEGAN v2.10)
    #                (2) 1/2 of SWRD is in 400-700nm band

    # compute Gmax limitation factors -----------------------------------------
    # ftemp -------------------------------------------------------------------
    BTexp   <- (p_specie$Tmax - p_specie$Topt) / (p_specie$Topt - p_specie$Tmin) 
    var_tmp <- (tem2loc - p_specie$Tmin) / (p_specie$Topt - p_specie$Tmin) * 
      ((p_specie$Tmax - tem2loc) / (p_specie$Tmax - p_specie$Topt))
    
    var_tmp[which(tem2loc < p_specie$Tmin)] <- 0 # ftemp = 0 if temperature < Tmin
    var_tmp[which(tem2loc > p_specie$Tmax)] <- 0 # ftemp = 0 if temperature > Tmax
    
    var_tmp <- var_tmp ** BTexp
    ftemp   <- apply(var_tmp, c(1,2), function(x,y) { max(x,y) }, y = p_specie$fmin) 

    # fvpd --------------------------------------------------------------------
    vapp   <- 611 * exp(17.502 * (tem2loc) / (tem2loc + 240.97))
    vpdkPa <- 1e-3 * vapp * (1 - srehloc)
    fvpd   <- p_specie$fmin + (1 - p_specie$fmin) * (p_specie$VPDmin - vpdkPa) / 
      (p_specie$VPDmin - p_specie$VPDmax)
    fvpd   <- apply(fvpd, c(1,2), function(x,y) { min(1, max(y,x)) }, y = p_specie$fmin)

    # flight ------------------------------------------------------------------
    flig <- 1 - exp(-p_specie$light_a * ppfd)

    # fsw ---------------------------------------------------------------------
    # Soil Moisture Index and Soil Water factor for PAW (crops)
    # all var (included PAWt) are alrdy between [0-1]
    SMI <- (swvl - PWP) / (FC - PWP)
    
    this_whichSMI <- which(SMI <= 0)
    if (length(this_whichSMI > 0)) { 
      fsmi[this_whichSMI] <- 0
    }
    
    this_whichSMI <- which((SMI >= 0) & (SMI < p_specie$PAWt))
    if (length(this_whichSMI) > 0) {
      fsmi[this_whichSMI] <- (SMI[this_whichSMI] / p_specie$PAWt)
    }
    
    this_whichSMI <- which(SMI >= p_specie$PAWt)
    if (length(this_whichSMI) > 0) {
      fsmi[this_whichSMI] <- 1
    }
    
    # Parameters from CHIMERE DEPO_SPEC:
    #charspec,rMx,rHx,rf0,rRwat
    #O3      48  0.01    1   2000
    factD <- sqrt(48 / 48)
    dHx <- 0.01 # from DEPO_SPEC for O3
    df0 <- 1    # from DEPO_SPEC for O3
    Rm <- 1e-2 / (dHx / 3000 + 100 * df0)

    # fsmi only used for crop species
    fswp2d <- array(1, dim = c(nlon,nlat)) # fsmi == 1 for tree (no limitation)
    if (specie %in% c("Wheat","Potato","Tomato")) {
      fswp2d <- fsmi 
    } # fswp
    
    mmax <- array(0, dim = c(nlon, nlat, 2))
    mmax[,,1] <- array(p_specie$fmin, dim = c(nlon,nlat)) # fmin2d
    mmax[,,2] <- ftemp * fvpd * fswp2d   
    
    mmin <- array(0, dim = c(nlon, nlat, 2))
    mmin[,,1] <- fphen[,,d]
    mmin[,,2] <- fO3 # fO3 of the previous day here 

    # deposition resistance & ozone uptake (ICP vegetation methodology --------
    O3Gsto_pody <- p_specie$gmax * apply(mmin, 1:2,min) * flig * apply(mmax, 1:2, max) / 41000
    Gsto_pody <- O3Gsto_pody
    Rc_pody <- 1 / (Gsto_pody + O3Gext)
    
    # compute O3 at z canopy using resistances Ra (background documents ICP VEGETATION)
    vkarm <- 0.41
    rough <- p_specie$zcanopy/10.
    Sc    <- 0.93 # Schmidt Number
    Pr    <- 0.71
    disp  <- p_specie$zcanopy * 2 / 3
    
    # downward condition on zcanopy (see Mapping Manual)
    if (p_specie$zcanopy < O3layerhght) {
      # init and compute similarity function 
      PhiKi1 <- PhiKi2 <- PhiKi3 <- array(0, dim = c(nlon,nlat))
      
      Ki1 <- (O3layerhght - disp) / obukloc
      Ki2 <- (p_specie$zcanopy - disp) / obukloc
      Ki3 <- rough / obukloc
      
      this_whichKi1Pls <- which(Ki1 >= 0)
      this_whichKi1Mns <- which(Ki1 < 0)
      
      this_whichKi2Pls <- which(Ki2 >= 0)
      this_whichKi2Mns <- which(Ki2 < 0)
      
      this_whichKi3Pls <- which(Ki3 >= 0)
      this_whichKi3Mns <- which(Ki3 < 0)
      
      # Similitary function for Heat PhiH(Ki)
      PhiKi1[this_whichKi1Pls] <- -5 * Ki1[this_whichKi1Pls]
      PhiKi1[this_whichKi1Mns] <- 2 * log((1 + ((1 - 16 * Ki1[this_whichKi1Mns]) 
                                                ** (1 / 4)) ** 2) / 2)
      
      PhiKi2[this_whichKi2Pls] <- -5 * Ki2[this_whichKi2Pls]
      PhiKi2[this_whichKi2Mns] <- 2 * log((1 + ((1 - 16 * Ki2[this_whichKi2Mns]) 
                                                ** (1 / 4)) ** 2) / 2)
      
      PhiKi3[this_whichKi3Pls] <- -5 * Ki3[this_whichKi3Pls]
      PhiKi3[this_whichKi3Mns] <- 2 * log((1 + ((1 - 16 * Ki3[this_whichKi3Mns]) 
                                                ** (1 / 4)) ** 2) / 2)
      
      RaZtZm <- 1 / (vkarm * ustaloc) * (log((O3layerhght - disp) / 
                                               (p_specie$zcanopy - disp)) - PhiKi1 + PhiKi2)
      
      RaDispZ0Z <- 1 / (vkarm * ustaloc) * (log((O3layerhght - disp) / 
                                                  (rough)) - PhiKi1 + PhiKi3)
      
      Rin <- 14 * dsai[d] * p_specie$zcanopy / ustaloc
      
      # Assumevertical wind profile close to neutral condition near the ground TEST
      RO3diff <- 2 / (vkarm * ustaloc) * (Sc / Pr) ** (2 / 3)  
     
      # Rsurf with lai
      Rsurf <- 1 / (dlai[d] * O3Gsto_pody + dsai[d] / Rext + 1 / (Rin + Rsoil))
      
      # downscalled ozone
      O3canopy <- O3hppb *  (1 - RaZtZm / (RaDispZ0Z + RO3diff + Rsurf))
      
    } else {
      # no correction of O3 level for specie for which zcanopy ~ = chimere first layer (trees)
      O3canopy <- O3hppb 
    }

    # Mapping Manual methodology for Rb
    rNu <- 0.15
    DH2O <- 0.25
    DH2O_Dx <- sqrt(48 / 18)
    factRb <- ((rNu / (DH2O * Pr)) * DH2O_Dx) ** (2 / 3) * 2 / vkarm
    Rb <- 1.3 * 150 * (p_specie$leafdimension * 1e-2 / 
                         (log(p_specie$zcanopy / rough) * ustaloc / vkarm)) ** 0.5
    
    # compute ozone uptake in fstoloclu ---------------------------------------
    # O3h in ppb, airmloc in molec/cm3 and fstoloclu in molec/cm3
    ratiogstolu <- O3Gsto_pody * Rc_pody / (Rb + Rc_pody)
    fstoloclu <- O3canopy * airmloc / 1e9 * ratiogstolu 

    #1e6 for cm3-> m3, 1e9 to get fstolu in nano mole ->fstolu (nmol)/m3
    fstolu <- fstoloclu / avogadro * 1e6 * 1e9 
    
    # compute fO3 -------------------------------------------------------------
    if (specie == "Wheat") { 
      #wheat inside => compute POD0
      fstoluPOD0 <- fstolu
      fstoluPOD0[which(fstolu < 0)] <- 0
      fstoluPOD0 <- fstoluPOD0 * 3600 / 1e6
      
      # select appropriate cells in the accumulation period
      this_whichFO3 <- which((swrd >= 50) & 
        ((cumtem > p_specie$fphen_1_ets) & (cumtem < p_specie$fphen_5_ets)))     
	  
      # compute additional pody0 for these cells at day d
      pody0[this_whichFO3] <- pody0[this_whichFO3] + fstoluPOD0[this_whichFO3] 
      # compute a new fO3 only for these cells, for other cells fO3 remains at old value.
      fO3[this_whichFO3] <- ((1 + (pody0[this_whichFO3] / p_specie$fO3cst) ** p_specie$fO3exp) ** (-1))
    }
  
    # Compute only for accumulation period: AOT0 for daylight 
    # (solar rad > 50 w/m2) select appropriate cells for the accumulation period
    this_whichDaylight <- which(swrd >= 50)    
    if (length(this_whichDaylight) > 0) {
      if (specie == "Potato") {
        # AOT0 is acumulated since Astart (Pleijel et al. 2007 &/or 2002
        this_whichFO3Potato <- which((swrd >= 50) & (((cumtem - PotatoCumtem) >= p_specie$fphen_1_ets) & 
          ((cumtem - PotatoCumtem) <= p_specie$fphen_2_ets)))
        
        AOT0[this_whichFO3Potato] <- AOT0[this_whichFO3Potato] + O3canopy[this_whichFO3Potato] / 1000  # convert ppb in ppm
        # compute a new fO3 only for these cells, for other cells fO3 remains at old value.
        fO3[this_whichFO3Potato] <- ((1 + (AOT0[this_whichFO3Potato] / p_specie$fO3AOT) ** p_specie$fO3exp) ** (-1))          
      }
    }

    # substraction of the Y threeshold and conversion of the flux into mmol m-2 PLA for one hour
    fstolu <- fstolu - y_pod
    fstolu[which(fstolu < 0)] <- 0
    fstolu <- fstolu * 3600 / 1e6
    
    # Computation of additional PODy at each hour for daylight (solar rad > 50 w/m2)
    this_whichNightlight <- which(swrd < 50)
    if (length(this_whichNightlight) > 0) {
      fstolu[this_whichNightlight] <- 0 # No accumulation during the night
      fsmi[this_whichNightlight]   <- 0 # No accumulation  during the night for the control output of smi
    }
    if (specie == "Wheat") {
      this_which <- which((cumtem > p_specie$fphen_1_ets) & (cumtem < p_specie$fphen_5_ets))
      if (length(this_which) > 0) {
        pody[this_which] <- pody[this_which] + fstolu[this_which]
      } # compute additional pody for these cells at day d
    }
    
    if (specie == "Potato") {
      this_which <- which(((cumtem - PotatoCumtem) >= p_specie$fphen_1_ets) & ((cumtem - PotatoCumtem) <= p_specie$fphen_2_ets))
      if (length(this_which) > 0) {
        pody[this_which] <- pody[this_which] + fstolu[this_which]
      } # compute additional pody for these cells at day d
    }
    
    if (specie == "Tomato") {
      this_which <- which((TomatoDegDay >= p_specie$Astart_ETS) & (TomatoDegDay <= p_specie$Aend_ETS))
      if (length(this_which) > 0) {
        pody[this_which] <- pody[this_which] + fstolu[this_which]
      }
    }
    
    if (specie %in% c("Beech", "TemperateOak")) {
      this_which <- which((Astart_FD <= d) & (Aend_FD >= d))
      if (length(this_which) > 0) {
        pody[this_which] <- pody[this_which] + fstolu[this_which]
      }
    }
    
    if (specie == "Spruce") {
      # select appropriate cells in the accumulation period for specie NorwaySpruce
      this_which <- which((tem2loc > p_specie$Tmin) & (tem2loc < p_specie$Tmax)) 
      if (length(this_which) > 0) {
        pody[this_which] <- pody[this_which] + fstolu[this_which]
      } # compute additional pody for these cells at day d
    }
    
     # no conditions for grass/forb parameterization
    if (specie == "Grass") {
      if ((d >= p_depo$depsgs) & (d <= p_depo$depegs)) {
        this_which <- which(O3canopy > 0)
        if (length(this_which) > 0) {
          pody[this_which] <- pody[this_which] + fstolu[this_which]
        }  
      }
    }
   
    if (length(this_which) > 0) {
      accperiod[this_which] = 1
    }

    if (heavy_output) {
      # write netcdf variable for each hour
      ncvar_put(ncout, def_pody, pody, start = c(1, 1, dh), count = c(-1, -1, 1))
      ncvar_put(ncout, def_O3, O3canopy, start = c(1, 1, dh), count = c(-1, -1, 1)) 
      ncvar_put(ncout, def_fO3, fO3, start = c(1, 1, dh), count = c(-1, -1, 1))
      ncvar_put(ncout, def_fsmi, fswp2d, start = c(1, 1, dh), count = c(-1, -1, 1))
      ncvar_put(ncout, def_ftemp, ftemp, start = c(1, 1, dh), count = c(-1, -1, 1))
      ncvar_put(ncout, def_flig, flig, start = c(1, 1, dh), count = c(-1, -1, 1))
      ncvar_put(ncout, def_fphen, fphen[,,d], start = c(1, 1, dh), count = c(-1, -1, 1))
      ncvar_put(ncout, def_fvpd, fvpd, start = c(1, 1, dh), count = c(-1, -1, 1))
      ncvar_put(ncout, def_accper, accperiod, start = c(1, 1, dh), count = c(-1, -1, 1))
    }
  }
} # loop day & hour

# write lon, lat and close hourly netcdf of pody & limitation functions
if (heavy_output) {
  ncvar_put(ncout, def_tim, as.character(times))
  ncvar_put(ncout, def_lon, lon)
  ncvar_put(ncout, def_lat, lat)
  nc_close(ncout)
}

# write small pody files 
save(lon, lat, pody, file = ofile)

# close input netcdf
nc_close(fic) 

