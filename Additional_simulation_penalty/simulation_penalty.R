options("rgdal_show_exportToProj4_warnings"="none") # suppress annoying warnings
library(raster)
library(rasterVis)
library(viridis)
library(ggplot2)
library(future)
library(furrr)
library(doRNG)
library(RStoolbox)
library(Rcpp)
library(geoR)
library(Rhabit)
library(TMB)
library(ctmm)
library(sp)
library(terra)
library(sf)
library(dplyr)
library(langevinSSM)
setwd("/Users/fannydupont/Desktop/PhD/LangevinSSM_Package/Additionnal_simulation_penalty")
## specify, compile, and load TMB model
compile("fitLangevin_maskpen_neg.cpp")
dyn.load(dynlib(paste0("fitLangevin_maskpen_neg")))
source("R/helper_functions.R")

simulate_boundary <- function(sca, n_vertices = sample(8:25, 1)) { #  sca controls the size of the polygon
  angles <- sort(runif(n_vertices, 0, 2*pi)) # random number of angles
  radius <- runif(n_vertices, 0.6, 1) * sca *2 # random angles
  x <- radius * cos(angles) # get cartesian
  y <- radius * sin(angles) # get cartesian
  poly <- rbind(cbind(x, y), c(x[1], y[1])) # last row is to close the polygon.
  st_sfc(st_polygon(list(poly)), crs = 3416)
}

model <- 1 # underdamped Langevin (model = 1) or overdamped Langevin (model = 0)
nsims <- 1 # number of simulations
nbAnimals <- 5 # number of tracks
obsPerAnimal <- 5000 # number of simulated locations per track
timeStep <- 0.01 # time scale of simulation
beta <- c(-4, 6, 5,-100) # resource selection coefficients
ncov <- length(beta) # number of spatial covariates
sigma <- 5 # speed parameter 
gamma <- 0.5 # friction parameter
psi <- 1 # error SD scaling parameter
scale_factor <-1
## sampling rate, missing data, and measurement error
samplingRate <- 1
propMissing <- 0
smaj.sd <- 3 # SD for semi-major error ellipse axis; smaj ~ abs(Normal(0,smaj.sd)): 0.05 (1%),0.5 (10%),1 (20%), 1.25 (25%),1.666667 (30%), 2 (40%), 2.5 (50%), 3 (60%)
smin.sd <- smaj.sd/2 # SD for semi-minor error ellipse axis; smin ~ abs(Normal(0,smin.sd))
eor <- c(0,180) # range for error ellipse orientation (in degrees from north); eor ~ Uniform(eor[1],eor[2])
measurementError <- list(smaj.sd=smaj.sd,smin.sd=smin.sd,eor=eor) # setting measurementError <- NULL adds no measurement error to observations


## specify scale, resolution, and spatial autocorrelation for covariates
sca <- 200
lim <- c(-1, 1, -1, 1)*sca
cropExtent <- extent(lim)
resol <- 1
covRange <- c(0.1,0.5)

## smooth gradient specifications
npoints <- 0
curweight <- 1/2
weights <- c(curweight,rep((1-curweight)/npoints,npoints)) # ignored if npoints=0
zetaScale <- 1
tracepar <- FALSE

langSim <- langTMB <- list()
covlist <- plot_true <- plot_error <- idx_best <- list()
gradient_at_true <- estpath <- list()
parMat<- matrix(NA,nrow=nsims,8)
if(model==0){
  colnames(parMat) <- c("sigma",paste0("beta",1:ncov),"BC_with","BC_without")
} else colnames(parMat) <- c("sigma",paste0("beta",1:ncov),"gamma","BC_with","BC_without")

set.seed(1998,kind="Mersenne-Twister",normal.kind = "Inversion")
prop_on_land <- array() # to count the proportion on location data (with error) on land
for(isim in 1:nsims){
  nllk = array()
  sim_init = list()
  parMat_init <- matrix(NA,nrow=2,8)
  if(model==0){
    colnames(parMat_init) <- c("sigma",paste0("beta",1:ncov),"BC_with","BC_without")
  } else colnames(parMat_init) <- c("sigma",paste0("beta",1:ncov),"gamma","BC_with","BC_without")
  
  cat("Simulation",isim,"\n")
  
  #######################
  ## Define covariates ##
  #######################
  message("   Generating covariates...")
  covlist[[isim]] <- list()
  
  for(i in 1:(ncov)) {
    irange <- runif(1,covRange[1],covRange[2])
     covlist[[isim]][[i]] <-  simCov(sca = sca, irange=irange, sigma2 = 0.1, kappa = 0.5)
     terra::crs(covlist[[isim]][[i]] ) <- "epsg:3416"
  }
 
  names(covlist[[isim]]) <- c("cov1","cov2","cov3","pen")
  template <- covlist[[isim]][[1]]
  terra::crs(template) <- "epsg:3416"
  
  # Start with all NA; only island cells will be filled with dist-to-interior,
  # ocean will be set to 0 at the end and boundary + dist-to-interior > 0.
  bathy <- template
  values(bathy) <- NA
  
  ext_center <- ext(template)
  center_x <- runif(1,-10,10)
  center_y <- runif(1,-10,10)
  
 # --- Island 1 (main, sca=40) --- Large island
  island1_sf <- simulate_boundary(sca = 30)                   # simulate a large random polygon (the "island") 
  island1_sf <- st_geometry(island1_sf)        
  island1_sf <- island1_sf + matrix(c(center_x, center_y))    # move the polygon to a specific center location: center_x, center_y
  island1_sf <- st_sf(geometry = island1_sf, crs = 3416)    
  island1_v  <- vect(island1_sf)                              # convert to terra SpatVector for raster operations
  
  island1_mask <- terra::rasterize(island1_v, template, field = 1)  # Values for the polygon: 1 inside, NA outside
  island1_mask[is.na(island1_mask)] <- 0       # replace NA (outside polygon) with 0 such that we get 1 = island, 0 = water
  
  inside_island1 <- ifel(island1_mask == 1, NA, 1)  # island = NA, water = 1
  # terra::distance needs NA where you want distances measured FROM.
  
  dist_from_edge1 <- terra::distance(inside_island1) # for each island cell, compute the distance to the nearest water cell (0 at shoreline, increases towards island center)

  
  # take the island's "distance to water" values and paste them into the bathymetry raster
  # cells outside the island are ignored (multiplying by island1_mask zeroes them out)
  
  dist_to_center1 <- dist_from_edge1    
  dist_to_center1[is.na(dist_to_center1)] <- 0  
  bathy[island1_mask == 1] <- values(dist_to_center1 * island1_mask)[island1_mask[] == 1] # positive values that increase from the shoreline inward, peaking at the center of the island

  
  # REPEAT THE SAME PROCEDURE FOR MANY ISLANDS ######
  # # --- Island 2  --- Large island
  other3_x <- runif(1,-40,-30)
  other3_y <- runif(1,180,190)

  island3_sf <- simulate_boundary(40)
  island3_sf <- st_geometry(island3_sf)
  island3_sf <- island3_sf + matrix(c(other3_x, other3_y))
  island3_sf <- st_sf(geometry = island3_sf, crs = 3416)
  island3_v <- vect(island3_sf)

  island3_mask <- terra::rasterize(island3_v, template, field = 1)
  island3_mask[is.na(island3_mask)] <- 0

  inside_island3 <- ifel(island3_mask == 1, NA, 1)
  dist_from_edge3 <- terra::distance(inside_island3)
  dist_to_center3 <- dist_from_edge3
  dist_to_center3[is.na(dist_to_center3)] <- 0
  bathy[island3_mask == 1] <- values(dist_to_center3 * island3_mask)[island3_mask[] == 1]

  # --- Small island A ---
  isl_a_x <- runif(1,-190,-180)
  isl_a_y <- runif(1,-190,-180)
  
  islandA_sf <- simulate_boundary(sca = 80)
  islandA_sf <- st_geometry(islandA_sf)
  islandA_sf <- islandA_sf + matrix(c(isl_a_x, isl_a_y))
  islandA_sf <- st_sf(geometry = islandA_sf, crs = 3416)
  islandA_v <- vect(islandA_sf)
  
  islandA_mask <- terra::rasterize(islandA_v, template, field = 1)
  islandA_mask[is.na(islandA_mask)] <- 0
  
  inside_islandA <- ifel(islandA_mask == 1, NA, 1)
  dist_from_edgeA <- terra::distance(inside_islandA)
  dist_to_centerA <- dist_from_edgeA
  dist_to_centerA[is.na(dist_to_centerA)] <- 0
  bathy[islandA_mask == 1] <- values(dist_to_centerA * islandA_mask)[islandA_mask[] == 1]
  
  # --- Small island B ---
  isl_b_x <- runif(1,-190,-180)
  isl_b_y <- runif(1,100,120)
  islandB_sf <- simulate_boundary(sca = 60)
  islandB_sf <- st_geometry(islandB_sf)
  islandB_sf <- islandB_sf + matrix(c(isl_b_x, isl_b_y))
  islandB_sf <- st_sf(geometry = islandB_sf, crs = 3416)
  islandB_v <- vect(islandB_sf)
  
  islandB_mask <- terra::rasterize(islandB_v, template, field = 1)
  islandB_mask[is.na(islandB_mask)] <- 0
  
  inside_islandB <- ifel(islandB_mask == 1, NA, 1)
  dist_from_edgeB <- terra::distance(inside_islandB)
  dist_to_centerB <- dist_from_edgeB
  dist_to_centerB[is.na(dist_to_centerB)] <- 0
  bathy[islandB_mask == 1] <- values(dist_to_centerB * islandB_mask)[islandB_mask[] == 1]
  
  # --- Small island C ---
  isl_c_x <- runif(1,100,120)
  isl_c_y <- runif(1,-100,-80)
  
  islandC_sf <- simulate_boundary(sca = 15)
  islandC_sf <- st_geometry(islandC_sf)
  islandC_sf <- islandC_sf + matrix(c(isl_c_x, isl_c_y))
  islandC_sf <- st_sf(geometry = islandC_sf, crs = 3416)
  islandC_v <- vect(islandC_sf)
  
  islandC_mask <- terra::rasterize(islandC_v, template, field = 1)
  islandC_mask[is.na(islandC_mask)] <- 0
  
  inside_islandC <- ifel(islandC_mask == 1, NA, 1)
  dist_from_edgeC <- terra::distance(inside_islandC)
  dist_to_centerC <- dist_from_edgeC
  dist_to_centerC[is.na(dist_to_centerC)] <- 0
  bathy[islandC_mask == 1] <- values(dist_to_centerC * islandC_mask)[islandC_mask[] == 1]
  
  # --- Small island D ---
  
  isl_d_x <- runif(1,160,190)
  isl_d_y <- runif(1,-100,-80)
  
  islandD_sf <- simulate_boundary(sca = 15)
  islandD_sf <- st_geometry(islandD_sf)
  islandD_sf <- islandD_sf + matrix(c(isl_d_x, isl_d_y))
  islandD_sf <- st_sf(geometry = islandD_sf, crs = 3416)
  islandD_v <- vect(islandD_sf)
  
  islandD_mask <- terra::rasterize(islandD_v, template, field = 1)
  islandD_mask[is.na(islandD_mask)] <- 0
  
  inside_islandD <- ifel(islandD_mask == 1, NA, 1)
  dist_from_edgeD <- terra::distance(inside_islandD)
  dist_to_centerD <- dist_from_edgeD
  dist_to_centerD[is.na(dist_to_centerD)] <- 0
  bathy[islandD_mask == 1] <- values(dist_to_centerD * islandD_mask)[islandD_mask[] == 1]
  
  # --- Small island E ---
  
  isl_e_x <- runif(1,110,130)
  isl_e_y <- runif(1,10,50)
  
  islandE_sf <- simulate_boundary(sca = 30)
  islandE_sf <- st_geometry(islandE_sf)
  islandE_sf <- islandE_sf + matrix(c(isl_e_x, isl_e_y))
  islandE_sf <- st_sf(geometry = islandE_sf, crs = 3416)
  islandE_v <- vect(islandE_sf)

  islandE_mask <- terra::rasterize(islandE_v, template, field = 1)
  islandE_mask[is.na(islandE_mask)] <- 0
  
  inside_islandE <- ifel(islandE_mask == 1, NA, 1)
  dist_from_edgeE <- terra::distance(inside_islandE)
  dist_to_centerE <- dist_from_edgeE
  dist_to_centerE[is.na(dist_to_centerE)] <- 0
  bathy[islandE_mask == 1] <- values(dist_to_centerE * islandE_mask)[islandE_mask[] == 1]
  
  # --- Small island G ---
  
  isl_g_x <- runif(1,100,120)
  isl_g_y <- runif(1,-150,-140)
  
  islandG_sf <- simulate_boundary(sca = 15)
  islandG_sf <- st_geometry(islandG_sf)
  islandG_sf <- islandG_sf + matrix(c(isl_g_x, isl_g_y))
  islandG_sf <- st_sf(geometry = islandG_sf, crs = 3416)
  islandG_v <- vect(islandG_sf)
  
  islandG_mask <- terra::rasterize(islandG_v, template, field = 1)
  islandG_mask[is.na(islandG_mask)] <- 0
  
  inside_islandG <- ifel(islandG_mask == 1, NA, 1)
  dist_from_edgeG <- terra::distance(inside_islandG)
  dist_to_centerG <- dist_from_edgeG
  dist_to_centerG[is.na(dist_to_centerG)] <- 0
  bathy[islandG_mask == 1] <- values(dist_to_centerG * islandG_mask)[islandG_mask[] == 1]

  # --- Small island H ---
  
  isl_h_x <- runif(1,0,15)
  isl_h_y <- runif(1,-150,-140)
  
  islandH_sf <- simulate_boundary(sca = 20)
  islandH_sf <- st_geometry(islandH_sf)
  islandH_sf <- islandH_sf + matrix(c(isl_h_x, isl_h_y))
  islandH_sf <- st_sf(geometry = islandH_sf, crs = 3416)
  islandH_v <- vect(islandH_sf)
  
  islandH_mask <- terra::rasterize(islandH_v, template, field = 1)
  islandH_mask[is.na(islandH_mask)] <- 0
  
  inside_islandH <- ifel(islandH_mask == 1, NA, 1)
  dist_from_edgeH <- terra::distance(inside_islandH)
  dist_to_centerH <- dist_from_edgeH
  dist_to_centerH[is.na(dist_to_centerH)] <- 0
  bathy[islandH_mask == 1] <- values(dist_to_centerH * islandH_mask)[islandH_mask[] == 1]
  
  # --- Small island I ---
  isl_i_x <- runif(1,160,170)
  isl_i_y <- runif(1,120,180)
  
  islandI_sf <- simulate_boundary(sca = 20)
  islandI_sf <- st_geometry(islandI_sf)
  islandI_sf <- islandI_sf + matrix(c(isl_i_x, isl_i_y))
  islandI_sf <- st_sf(geometry = islandI_sf, crs = 3416)
  islandI_v <- vect(islandI_sf)
  
  islandI_mask <- terra::rasterize(islandI_v, template, field = 1)
  islandI_mask[is.na(islandI_mask)] <- 0
  
  inside_islandI <- ifel(islandI_mask == 1, NA, 1)
  dist_from_edgeI <- terra::distance(inside_islandI)
  dist_to_centerI <- dist_from_edgeI
  dist_to_centerI[is.na(dist_to_centerI)] <- 0
  bathy[islandI_mask == 1] <- values(dist_to_centerI * islandI_mask)[islandI_mask[] == 1]
  
  # --- Small island J ---
  
  isl_j_x <- runif(1,-120,-100)
  isl_j_y <- runif(1,-30,-20)
  
  islandJ_sf <- simulate_boundary(sca = 20)
  islandJ_sf <- st_geometry(islandJ_sf)
  islandJ_sf <- islandJ_sf + matrix(c(isl_j_x, isl_j_y))
  islandJ_sf <- st_sf(geometry = islandJ_sf, crs = 3416)
  islandJ_v <- vect(islandJ_sf)
  
  islandJ_mask <- terra::rasterize(islandJ_v, template, field = 1)
  islandJ_mask[is.na(islandJ_mask)] <- 0
  
  inside_islandJ <- ifel(islandJ_mask == 1, NA, 1)
  dist_from_edgeJ <- terra::distance(inside_islandJ)
  dist_to_centerJ <- dist_from_edgeJ
  dist_to_centerJ[is.na(dist_to_centerJ)] <- 0
  bathy[islandJ_mask == 1] <- values(dist_to_centerJ * islandJ_mask)[islandJ_mask[] == 1]
  
  # --- Small island K ---
  isl_k_x <- runif(1,-10,10)
  isl_k_y <- runif(1,85,95)
  
  islandK_sf <- simulate_boundary(sca = 15)
  islandK_sf <- st_geometry(islandK_sf)
  islandK_sf <- islandK_sf + matrix(c(isl_k_x, isl_k_y))
  islandK_sf <- st_sf(geometry = islandK_sf, crs = 3416)
  islandK_v <- vect(islandK_sf)
  
  islandK_mask <- terra::rasterize(islandK_v, template, field = 1)
  islandK_mask[is.na(islandK_mask)] <- 0
  
  inside_islandK <- ifel(islandK_mask == 1, NA, 1)
  dist_from_edgeK <- terra::distance(inside_islandK)
  dist_to_centerK <- dist_from_edgeK
  dist_to_centerK[is.na(dist_to_centerK)] <- 0
  bathy[islandK_mask == 1] <- values(dist_to_centerK * islandK_mask)[islandK_mask[] == 1]
  
  # --- Small island L ---
  isl_l_x <- runif(1,-70,-60)
  isl_l_y <- runif(1,70,80)
  
  islandL_sf <- simulate_boundary(sca = 10)
  islandL_sf <- st_geometry(islandL_sf)
  islandL_sf <- islandL_sf + matrix(c(isl_l_x, isl_l_y))
  islandL_sf <- st_sf(geometry = islandL_sf, crs = 3416)
  islandL_v <- vect(islandL_sf)
  
  islandL_mask <- terra::rasterize(islandL_v, template, field = 1)
  islandL_mask[is.na(islandL_mask)] <- 0
  
  inside_islandL <- ifel(islandL_mask == 1, NA, 1)
  dist_from_edgeL <- terra::distance(inside_islandL)
  dist_to_centerL <- dist_from_edgeL
  dist_to_centerL[is.na(dist_to_centerL)] <- 0
  bathy[islandL_mask == 1] <- values(dist_to_centerL * islandL_mask)[islandL_mask[] == 1]
  
  # --- Small island M ---
  isl_m_x <- runif(1,-80,-70)
  isl_m_y <- runif(1,40,45)
  
  islandM_sf <- simulate_boundary(sca = 10)
  islandM_sf <- st_geometry(islandM_sf)
  islandM_sf <- islandM_sf + matrix(c(isl_m_x, isl_m_y))
  islandM_sf <- st_sf(geometry = islandM_sf, crs = 3416)
  islandM_v <- vect(islandM_sf)
  
  islandM_mask <- terra::rasterize(islandM_v, template, field = 1)
  islandM_mask[is.na(islandM_mask)] <- 0
  
  inside_islandM <- ifel(islandM_mask == 1, NA, 1)
  dist_from_edgeM <- terra::distance(inside_islandM)
  dist_to_centerM <- dist_from_edgeM
  dist_to_centerM[is.na(dist_to_centerM)] <- 0
  bathy[islandM_mask == 1] <- values(dist_to_centerM * islandM_mask)[islandM_mask[] == 1]
  
  # Normalise all island cells by global max dist-to-interior.
  # ALL remaining NA cells — ocean AND raster boundary (x/y = +/-200) — become 0 (dark purple).
  global_max_dist <- max(
    max(values(dist_to_center1), na.rm = TRUE),
    max(values(dist_to_centerA), na.rm = TRUE),
    max(values(dist_to_centerB), na.rm = TRUE),
    max(values(dist_to_centerC), na.rm = TRUE),
    max(values(dist_to_centerD), na.rm = TRUE),
    max(values(dist_to_centerE), na.rm = TRUE),
    max(values(dist_to_centerG), na.rm = TRUE),
    max(values(dist_to_centerH), na.rm = TRUE),
    max(values(dist_to_centerI), na.rm = TRUE),
    max(values(dist_to_centerJ), na.rm = TRUE),
    max(values(dist_to_centerK), na.rm = TRUE),
    max(values(dist_to_centerL), na.rm = TRUE),
    max(values(dist_to_centerM), na.rm = TRUE)
    
  )
  bathy_vals <- values(bathy)
  bathy_vals[!is.na(bathy_vals)] <- bathy_vals[!is.na(bathy_vals)] / global_max_dist
  bathy_vals[is.na(bathy_vals)] <- 0  # ocean cells get 0 (dark purple, no penalty)
  
  # Set raster boundary cells (x/y = +/-200) to 1 — same as deepest island interior,
  # so beta[4]*1 = -100 strongly repels animals from the map edge
  # Set a 5-unit buffer around the raster edge (x/y from ±190 to ±200) to 1
  xy <- xyFromCell(bathy, 1:ncell(bathy))
  boundary_mask <- xy[,1] <= (-sca + 20) | xy[,1] >= (sca - 20) |
    xy[,2] <= (-sca + 10) | xy[,2] >= (sca - 20)
  bathy_vals[boundary_mask] <- 0.9
  
  values(bathy) <- bathy_vals
  
  vals_mat <- matrix(terra::values(bathy), nrow = nrow(bathy), ncol = ncol(bathy), byrow = FALSE)
  xgrid <- seq(xmin(bathy) + 0.5, xmax(bathy) - 0.5, length.out = ncol(bathy))
  ygrid <- seq(ymin(bathy) + 0.5, ymax(bathy) - 0.5, length.out = nrow(bathy))
  
  bathy_rhabit <- list(x = xgrid, y = ygrid, z = vals_mat)
  covlist[[isim]][[4]] <- rast(bathy_rhabit)
  terra::crs(covlist[[isim]][[4]]) <- "epsg:3416"
  
  
  # Compute utilization distribution ######
  UD_with <- langevinSSM::getUD(spatialCovs=covlist[[isim]], beta=c(-4, 6, 5,-100))
  UD_without <-  langevinSSM::getUD(spatialCovs=covlist[[isim]], beta= c(-4, 6, 5,0))

  plot(covlist[[isim]][[4]])
  plot(UD_with)
  
  
  terra::crs(covlist[[isim]][[4]]) <- "epsg:3416"
  names(covlist[[isim]]) <- c("cov1","cov2","cov3","cov4")

  spatialCovs <-  covlist[[isim]]
  
  names(spatialCovs) <- c("cov1","cov2","cov3","cov4")
  
  # Combine list of SpatRasters into one multi-layer SpatRaster
  raster_stack <- terra::rast(spatialCovs)
  
  raster_data <- list(
    raster_vals    = array(as.numeric(terra::values(raster_stack)),
                           dim = c(terra::nlyr(raster_stack),
                                   nrow(raster_stack),
                                   ncol(raster_stack))),
    raster_coords     = terra::crds(raster_stack),
    raster_resolution = terra::res(raster_stack),
    raster_extent     = as.vector(terra::ext(raster_stack)),
    n_covs            = terra::nlyr(raster_stack)
  )
  
  
  # simulate "high resolution" tracks
  message("   Simulating tracks...")
  par <- list(beta=c(-4, 6, 5,-100),sigma=sigma,gamma=gamma,psi=psi)
  
  langSim[[isim]] <- tryCatch(stop(),error=function(e) e)
  attempts <- 0
  
  while(inherits(langSim[[isim]],"error") && attempts < 200){
    attempts <- attempts + 1
    
    langSim[[isim]] <- tryCatch(langevinSSM::simLangevin(
      model = "underdamped",
      nbAnimals = nbAnimals,
      obsPerAnimal = obsPerAnimal,
      timeStep = timeStep,
      par=par,
      spatialCovs = spatialCovs,
      measurementError = measurementError
    ),error=function(e) e)
    if(inherits(langSim[[isim]],"error")){
      message("    Retrying Simulation ",isim,": ",langSim[[isim]]$message)
    }
  } 
  if(inherits(langSim[[isim]],"error")){
    message("    Skipping Simulation ",isim," after 200 failed attempts.")
    next
  }
  
  
  # Plot tracks with boundary covariate
  track_df <- data.frame(x = langSim[[isim]]$mu.x, y = langSim[[isim]]$mu.y)

  
  bathy_vals <- terra::extract(bathy, track_df[, c("x", "y")])
  prop_on_land[isim] <-   mean(bathy_vals[, 2] > 0, na.rm = TRUE)
  
  # Convert SpatRaster to data frame for ggplot
  UD_df <- as.data.frame(UD_with, xy = TRUE)
  colnames(UD_df) <- c("x", "y", "value")
  
  track_df <- data.frame(x = langSim[[isim]]$mu.x, 
                         y = langSim[[isim]]$mu.y)
  
  plot_true[[isim]] <- ggplot() +
    geom_raster(data = UD_df, aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_c(name = expression(c[4])) +
    geom_point(data = track_df, aes(x = x, y = y), 
               color = "red", size = 1) +
    coord_equal() +
    theme_minimal()
  
  
  track_df <- data.frame(x = langSim[[isim]]$x, y = langSim[[isim]]$y)
  plot_error[[isim]] <- ggplot() +
    geom_raster(data = UD_df, aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_c(name = expression(c[4])) +
    geom_point(data = track_df, aes(x = x, y = y), 
               color = "red", size = 1) +
    coord_equal() +
    theme_minimal()
  
  # See code comments in Case_Study_code.Rmd for more explanation
  raster_stack <- covlist[[isim]]$cov4 # d2_water
  
  mask_data <- list(
    mask_vals    = array(as.numeric(terra::values(raster_stack)),
                         dim = c(terra::nlyr(raster_stack),
                                 nrow(raster_stack),
                                 ncol(raster_stack))),
    mask_coords     = terra::crds(raster_stack),
    mask_resolution = terra::res(raster_stack),
    mask_extent     = as.vector(terra::ext(raster_stack))
  )
  
  
  # subsample data
  probs <- rep(1,nrow(langSim[[isim]]))
  probs[cumsum(c(1,table(langSim[[isim]]$ID)[1:(nbAnimals-1)]))] <- 1.e+10
  subDat <- subSampleData(langSim[[isim]], samplingRate = samplingRate, propMissing = propMissing)
  
  data <- list(model=model,
               Y=t(subDat[,c("mu.x","mu.y")])/scale_factor,
               dt=subDat$dt)
  
  probs <- rep(1,nrow(subDat))
  probs[cumsum(c(1,table(subDat$ID)[1:(nbAnimals-1)]))] <- 0
  data$Y[,sample.int(nrow(subDat),nrow(subDat)*propMissing,prob=probs,replace=FALSE)] <- NA
  
  data$isd <- as.numeric(!is.na(data$Y[1,]))
  data$obs_mod <- rep(NA,ncol(data$Y))
  data$obs_mod[data$isd==1] <- 1
  data$ID <- subDat$id
  data$time <- subDat$date
  data$nbStates <- 1
  data$nbObs <- rep(1,ncol(data$Y))
  data$smaj <- subDat$smaj / scale_factor
  data$smin <- subDat$smin / scale_factor
  data$eor <- subDat$eor
  data$K <- matrix(NA,ncol(data$Y),2)
  data$scale_factor <- scale_factor
  data <- c(data,raster_data,mask_data)
  data$smoothGradient <- ifelse(npoints>0,1,0)
  data$weights <- weights
  data$zetaScale <- zetaScale
  
  notNA <- which(!is.na(data$Y[1,]))


langTMB_rep <- list()
for(j in (1:2)){
  if(j ==1){
    par <- initialValues(subDat,model="underdamped",spatialCovs=spatialCovs)
    parm <- list(
                 beta = matrix(c(par$beta[1:3],-1e4),1,4),
                 log_gamma = log(par$gamma),
                 log_sigma = log(par$sigma),
                 vel = t(par$vel),
                 mu = t(par$mu),
                 l_delta = 0,
                 l_gamma = matrix(0,1,2),
                 l_psi = log(psi),
                 l_tau = c(0,0),
                 l_rho_o = 0)
}else{
  sigma <- 5       # movement scale
  gamma <- 0.5       # friction/velocity persistence
  betapar <- c(-4, 6, 5,log(100))           # resource selection coefficients
  psi <- 1                         # measurement error SD scaling
  
  init.mu <- as.matrix(rbind(data$Y[1,], data$Y[2,]))
  # Compute initial velocities for each animal
  init.v_mu <- matrix(0, 2, ncol(data$Y))
  for(i in 1:nbAnimals){
    aInd <- which(data$id == unique(data$id)[i])
    init.v_mu[1,aInd[-1]] <- diff(init.mu[1,aInd]) / 
      (data$dt[aInd[-1]] * exp(-gamma * data$dt[aInd[-1]]))
    init.v_mu[2,aInd[-1]] <- diff(init.mu[2,aInd]) / 
      (data$dt[aInd[-1]] * exp(-gamma * data$dt[aInd[-1]]))
  }
  
  
  # Parameter list for TMB
  parm <- list(
    log_sigma = log(sigma),
    beta = matrix(betapar,1,length(betapar)),
    mu = init.mu,           # initial locations
    vel = init.v_mu,       # initial velocities
    log_gamma = log(gamma),
    l_delta = 0,
    l_gamma = matrix(0,1,2),
    l_psi = log(psi),
    l_tau = c(0,0),
    l_rho_o = 0
  )
}

re <- c("mu", "vel")        # treat locations and velocities as random effects
map_inner <- lapply(parm[names(parm)[!names(parm) %in% c("mu","vel")]], function(x) factor(rep(NA,length(x))))

obj1 <- TMB::MakeADFun(
  data,
  parm,
  map = map_inner,
  random = re,
  DLL = "fitLangevin_maskpen_neg",
  hessian = TRUE,
  silent = TRUE,
  inner.control = list(maxit = 10000, trace=0)
)

# obj1$fn()

smoothed_pars <- obj1$env$parList()

parm$mu <- smoothed_pars$mu
parm$vel <- smoothed_pars$vel

map <- list(
  l_delta  = factor(NA),
  l_gamma  = factor(c(NA, NA)),
  l_rho_o  = factor(NA),
  l_tau    = factor(c(NA, NA)),
  l_psi = factor(NA)
)
# Create AD function for TMB
message("   Fitting model...")
obj2 <- MakeADFun(
  data,
  parm,
  map = map,
  random = re,
  DLL = "fitLangevin_maskpen_neg",
  hessian = TRUE,
  silent = TRUE,
  inner.control = list(maxit = 10000, trace=0)
)

obj2$env$tracepar <- tracepar


langTMB_rep[[j]] <- tryCatch(
      do.call("nlminb", args = list(
        start = obj2$par,
        objective = obj2$fn,
        gradient = obj2$gr,
        control=list(trace=10,iter.max=10000,eval.max=10000))),
      error   = function(e) { message("Error in iteration ",   j, ": ", e$message); return(NA) },
      warning = function(w) { message("Warning in iteration ", j, ": ", w$message); return(NA) }
)
    
if(sum(is.na(langTMB_rep[[j]]))==0){
    message("   Calculating SEs...")
    rep <- obj2$report()
    sdrep <- sdreport(obj2)
    langTMB_rep[[j]]$sdreport <- summary(sdrep,"report")
      
    nllk[j] <- ifelse(langTMB_rep[[j]]$convergence==1, NA, obj2$fn())
    parMat_init[j,1:6] <- langTMB_rep[[j]]$par
    parMat_init[j,"sigma"] <- rep$sigma
    if(model==1) parMat_init[j,"gamma"] <- rep$gamma
      
    estUD <- langTMB_rep[[j]]$par[2]*spatialCovs$cov1 +
      langTMB_rep[[j]]$par[3]*spatialCovs$cov2 +
      langTMB_rep[[j]]$par[4]*spatialCovs$cov3 #+
      #langTMB_rep[[j+1]]$par[6]*spatialCovs$cov4
      
      estUD_exp      <- exp(values(estUD))
      UD_with_exp    <- exp(values(UD_with))
      UD_without_exp <- exp(values(UD_without))
      
      estUD_exp      <- estUD_exp / sum(estUD_exp)
      UD_with_exp    <- UD_with_exp / sum(UD_with_exp)
      parMat_init[j,"BC_with"] <- sum(sqrt(estUD_exp * UD_with_exp))
      
      UD_without_exp <- UD_without_exp / sum(UD_without_exp)
      parMat_init[j,"BC_without"] <- sum(sqrt(estUD_exp * UD_without_exp))
      print(parMat_init[j,"BC_without"])
    } else {
      nllk[j] <- NA
    }
  }
  
  idx = which.min(nllk)
  if(length(idx)>0){
    parMat[isim,] <- parMat_init[idx,]
    langTMB[[isim]] <- langTMB_rep[[idx]]
    idx_best[[isim]] <- idx
  } else {
    parMat[isim,] <- rep(NA,ncol(parMat))
    langTMB[[isim]] <- NULL
  }
}
# 
# save(parMat,beta,sigma,gamma,prop_on_land,obsPerAnimal,langSim,langTMB,covlist,psi,timeStep,
#      samplingRate,propMissing,smaj.sd,scale_factor,resol,sca,npoints,covRange,
#      file=paste0("Results_penalty_new/",ifelse(model==1,"underdamped","overdamped"),
#                  "Sim_BC_sigmapen_withcov_inflate_",1,
#                  "_nbAnimals_",nbAnimals,
#                  "_obsPerAnimal_",obsPerAnimal,
#                  "_timeStep_",timeStep,
#                  "_sigma_",sigma,
#                  "_gamma_",gamma,
#                  "_beta_",paste0(beta,collapse="_"),
#                  "_samplingRate_",samplingRate,
#                  "_propMissing_",propMissing,
#                  "_M_",smaj.sd,
#                  "_npoints_",npoints,
#                  "_psi_",psi,
#                  ".RData"))

