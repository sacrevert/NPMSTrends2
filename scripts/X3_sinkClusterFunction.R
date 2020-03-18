#######################################################################################
## Function for preparing species datasets for the JAGS model and then running model ##
#######################################################################################
runModels_v5c_CLUS <- function(i, dat, file = file, bHab = bHab,
                               n.chains = 6, n.adapt = 100, n.iter = 5000, thin = 5){ #### PROBABLY NEED TO ALTER FILE LOCATION FOR CLUSTER ####
  x <- dat[[1]]
  #x <- (x %>% group_by(plot_id) %>% filter(any(!is.na(dominUnify))) %>% as.data.frame()) # added in v1: keep plots with at least one non NA only
  tryCatch( {
    if ( is.factor(x$date.x) ) { x$date.x <- as.Date(as.character(x$date.x), format = "%d/%m/%Y") }
    if ( is.character(x$grazing) ) {x$grazing[is.na(x$grazing)] <- "Absent" }
    if ( is.factor(x$grazing) ) {x$grazing <- as.character(x$grazing); x$grazing[is.na(x$grazing)] <- "Absent"}
    x$grazing <- factor(x$grazing, ordered = T, levels = c("Absent", "Low", "Moderate", "High"))
    #levels(x$cutting) <- "cut/mow" # not needed yet!
    x$title <- factor(x$title, ordered = T, levels = c("Wildflower survey", "Indicator survey", "Inventory survey"))
    x$year <- format(x$date.x, "%Y")
    x <- x[order(x$year, x$plot_id), ]
    uniPlots <- unique(x$plot_id) # unique plot IDs - 31/12/2018 = 933
    # create unique plot index (useful for cross-referencing later on?)
    plotIndex <- data.frame(plot = uniPlots, index = 1:length(uniPlots))
    N <- length(uniPlots)
    # check that all years are in the range of 2015:System time
    if ( max(unique(format(x$date.x, "%Y"))) > format(Sys.time(), "%Y") ) {
      print(paste("Error. Year violation: ", names(dat)))
    } else {
      print(paste("Years ok: ", names(dat)))
    }
    Y <- length(unique(format(x$date.x, "%Y")))
    # All Domin data, including zeros
    yOrig <- x$dominUnify
    x2 <- x$grazing
    x3 <- x$title
    n.Plot.pos <- length(x$dominUnify[x$dominUnify !='0' & !is.na(x$dominUnify)]) # 300
    cpos.Cens <- rep(1, n.Plot.pos)
    cpos.Latent <- rep(NA, n.Plot.pos)
    t <- c(1e-4, 0.01,
           0.01, 0.03,
           0.03, 0.05,
           0.05, 0.1,
           0.1, 0.25,
           0.25, 0.33,
           0.33, 0.5,
           0.5, 0.75,
           0.75, 0.95, 
           0.95, 0.9999,
           1e-4,0.9999) #this row is for presences with unknown cover (given the value 11 in the Domin scale)
    tdf <- as.data.frame(matrix(t, nrow = 11, ncol = 2, byrow = TRUE))
    colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
    tdf$int <- c(1,2,3,4,5,6,7,8,9,10,11)
    spPos <- x[x$dominUnify !='0' & !is.na(x$dominUnify),]
    # add row number to allow resorting after merge (merge with sort = F does not actually do what we want, could also used plyer::join)
    spPos$indexPos <- 1:nrow(spPos) 
    spPos <- merge(spPos, tdf, by.x = "dominUnify", by.y = "int", all.x = T, all.y = F)
    spPos <- spPos[order(spPos$indexPos),]
    lims <- spPos[,c("L","U")] # has the lower/upper cutpoints for all intervals
    spPos <- merge(spPos, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
    spPos <- spPos[order(spPos$indexPos),]
    plot <- spPos$index
    spPos$year <- as.factor(spPos$year)
    levels(spPos$year) <- 1:length(unique(format(x$date.x, "%Y")))
    year <- spPos$year # an indicator linking a percentage cover observation to its parent year
    V2 <- nrow(x) # the total number of visits (samples), irrespective of whether there is a positive cover for a species or not
    #plotZ <- match(Achi_mill_PAN$plot_id, uniPlots)  # an indicator linking a visit to its parent (spatially unique) plot
    x <- merge(x, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
    x <- x[order(x$year, x$plot_id), ]
    plotZ <- x$index
    x$year <- as.factor(x$year)
    levels(x$year) <- 1:length(unique(format(x$date.x, "%Y")))
    yearZ <- x$year # an indicator linking a visit to its parent year
    y <- y1 <-  numeric()
    for (j in 1:nrow(x)) {
      if(is.na(x$dominUnify[j])){ 
        y <- NA
      } else if (x$dominUnify[j]==0) { 
        y <- 0
      } else {
        y <- 1
      }
      y1 <- c(y1,y)
    }
    y <- y1 # visit-level detection history (binary)
    # Prepare data for JAGS
    yOrig[yOrig==11] <- NA # remove these for the purposes of using yOrig as a covar in detection model
    yOrig[is.na(yOrig)] <- 0 # no cover = zero for covar purposes
    y.ind <- yOrig # also use Orig to create indicator to indicate whether there is a cover value or one needs creating
    y.ind[y.ind >= 1] <- 1 # yOrig cover value present
    y.ind[is.na(y.ind)] <- 0 # yOrig cover value present
    Data <- list(N = N,
                 Y = Y,
                 n.Plot.pos = n.Plot.pos,
                 cposCens = cpos.Cens, # indicator (is censored?)
                 cposLatent = cpos.Latent, # NA values for latent observations
                 lims = lims,
                 cuts = c(1e-16,1e-4,0.01,0.03,0.05,0.1,0.25,0.33,0.5,0.75,0.95,0.9999),
                 year = year,
                 V2 = V2,
                 plotZ = plotZ,
                 yearZ = yearZ,
                 y = y,
                 x2 = as.numeric(x2), # grazing
                 x3 = as.numeric(x3), # survey level
                 yOrig = yOrig, # make sure yOrig = 11 is deleted (used as covar on detection)
                 y.ind = y.ind,
                 g2Levs = length(unique(x2)), # sometimes some levels might be absent
                 g3Levs = length(unique(x3))) # sometimes some levels might be absent
    zinit <- matrix(1, nrow = N, ncol = Y)
    inits.fn <- function() list(z = zinit,
                                sdC = runif(1,0,10),
                                sd0 = runif(1,0,10),
                                sd.m = runif(1,0,10),
                                sd.g2 = runif(1,0,10),
                                sd.g3 = runif(1,0,10),
                                mean.m = runif(1,0,1),
                                mean.p = runif(1,0,1),
                                mean.g2 = runif(1,0,1),
                                mean.g3 = runif(1,0,1),
                                phi = rep(runif(1,1,10), Y),
                                y.hat = rep(runif(1,0,1), V2),
                                g1 = runif(1,0,3),
                                # note extra initial value for dominUnify = 11 (represents a presence with unknown cover)
                                cposLatent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify]
    )
    #for ref only
    cPos.Init <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify]
    ### MAKE SURE YOU HAVE GIVEN THE RIGHT MODEL SCRIPT TO THE FUNCTION ###
    log <- file(paste("scripts/error_file_", i, "_", bHab, ".log", sep = ""), open = "wt")
    sink(log, append = TRUE, type = "message")
    jagsModel <- NULL
    attempt <- 0
    while( is.null(jagsModel) && attempt <= 3 ) {
      attempt <- attempt + 1
      try( {jagsModel <- rjags::jags.model(file = file, data = Data, 
                                       inits = inits.fn, n.chains = n.chains, n.adapt = n.adapt)}, 
                     silent = F)
      }
    ### MAKE SURE YOU HAVE THE RIGHT MODEL SCRIPT ###
    # Specify parameters for which posterior samples are saved
    #para.names <- c('psi')
    para.names <- c('mC', 'mPsi', 'mu', 'annOcc', 'avgOcc')
    # Continue the MCMC runs with
    attempt2 <- 0
    while(attempt2 <= 3){
      attempt2 <- attempt2 + 1
      try( { update(jagsModel, n.iter = n.iter)}, # burn in!
        silent = F)
    }
    attempt3 <- 0
    samples <- NULL
    while( is.null(samples) && attempt3 <= 3 ) {
      attempt3 <- attempt3 + 1
      try ( {samples <- rjags::coda.samples(jagsModel, variable.names = para.names, n.iter = n.iter, thin = thin)},
            silent = F)
    }
    sink(type = "message") ## sink any captured error messages to file
    ## Inspect results
    out <- summary(samples)
    mu_sd <- out$stat[,1:2] #make columns for mean and sd
    q <- out$quantile[,c(1,3,5)] #make columns for median and CI
    tableOut <- as.data.frame(cbind(mu_sd,q)) #make table
    ##################################
    ## add DIC or similar as an output
    ##################################
    return(list(x,tableOut,samples)) },
    error = function(err) print(x)) # print x, which is species name only if there are no samples
}