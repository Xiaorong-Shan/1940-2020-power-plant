library( data.table)
library( raster)
# library( spBayes)
library( disperseR)
library( ggplot2)
library( viridis)
library( lubridate)
library( mgcv)
# library( xgboost)
library( pbmcapply)

`%ni%` <- Negate(`%in%`) 

#======================================================================#
## ddm to grid raster
#======================================================================#
ddm_to_zip <- function( ddm_coal_file,
                        Year,
                        avg.period = 'year'){
  
  p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
  
  #read data, 
  ddm_coal <- fread(ddm_coal_file)
  
  # melt, extract year, month, day
  ddm_coal.m <- melt( ddm_coal, id.vars = c( 'X', 'Y'),
                      variable.name = 'date.in', value.name = 'coal_pm25')
  ddm_coal.m[, `:=` ( date.in = as.Date( date.in, format = '%m/%d/%y'))]
  ddm_coal.m[, `:=` ( year.in = year( date.in),
                      month.in = month( date.in))]
  
  # remove blow up values (greater than 30)
  ddm_coal.m[ coal_pm25 < 0 | coal_pm25 > 30, coal_pm25 := NA]
  
  #rasterize as brick
  names.ddm <- unique( ddm_coal.m$date.in)
  ddm_coal.b <- brick( lapply( names.ddm,
                               function( name, dt.m){
                                 r <- rasterFromXYZ(dt.m[ date.in == name,
                                                          .(x = X, y = Y, z = coal_pm25)],
                                                    crs = CRS(p4s))
                                 names( r) <- name
                                 return( r)
                               }, ddm_coal.m))
  
  # fill NA's with linear interpolation across days
  ddm_coal.b <- approxNA( ddm_coal.b, rule=2)
  names( ddm_coal.b) <- names.ddm
  
  # take monthly averages
  if( avg.period == 'month'){
    ddm_coal.month <- lapply( 1:12,
                              function( m, ddm_raster.b){
                                names.dates <- as.Date( gsub( '\\.', '-', gsub( 'X', '', names( ddm_raster.b))))
                                id <- which( month( names.dates) == m)
                                
                                ddm_coal.mon <- mean( subset( ddm_raster.b, id))
                                
                                return( ddm_coal.mon)
                              }, ddm_coal.b)
    
    ddm_coal.z <- brick( ddm_coal.month)
    names( ddm_coal.z) <- paste( Year, formatC( 1:12, width = 2, flag = '0'), sep = '.')
    
    # take annual averages
  } else if( avg.period == 'year'){
    # take annual average
    ddm_coal.z <- mean( ddm_coal.b)
    names( ddm_coal.z) <- Year
    
  }
  
  return( ddm_coal.z)
}


#======================================================================#
## functions to get meteorology data
# download the necessary met files, 20th century reanalysis
#======================================================================#
downloader.fn <- function( filename,
                           destination = file.path('~', 'Dropbox', 'Harvard', 'RFMeval_Local', 
                                                   'Comparisons_Intermodel', 'Global_meteorology'),
                           dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR')){
  if( length( dataset) > 1)
    dataset <- dataset[1]
  fileloc <- file.path( destination, dataset)
  
  
  # create directory to store in
  dir.create( fileloc, 
              recursive = T, 
              showWarnings = F)
  
  # name variable, filenames
  varname_NOAA <- gsub( "\\..*", "", filename)
  file_NOAA <- file.path( fileloc, filename)
  
  # define URL
  if( dataset == '20thC_ReanV2c'){
    # https://www.esrl.noaa.gov/psd/data/gridded/data.20thC_ReanV2c.monolevel.mm.html
    url_NOAA <- paste0( "ftp://ftp.cdc.noaa.gov/Datasets/20thC_ReanV2c/Monthlies/gaussian/monolevel/", filename)
  } else if( dataset == 'ncep.reanalysis.derived'){
    url_NOAA <- paste0( "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/", filename)
  }else if( dataset == 'NARR'){
    url_NOAA <- paste0( "ftp://ftp2.psl.noaa.gov/Datasets/NARR/Monthlies/monolevel/", filename) 
		      #  "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Monthlies/monolevel/", filename)
  }
  
  if( !file.exists( file_NOAA))
    download.file( url = url_NOAA,
                   destfile = file_NOAA)
  
  hpbl_rasterin <- brick( x = file_NOAA, 
                          varname = varname_NOAA)
  
  return( hpbl_rasterin)
  
}

#======================================================================#
## functions to get meteorology data
# extract the year of interest, average by year or return months
#======================================================================#
extract_year.fn <- function( raster.in = list.met[[1]],
                             year.in = 2005,
                             avg.period = 'year',
                             dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR')){
  
  # default to 20th cent reanalysis
  if( length( dataset) > 1){
    dataset <- dataset[1]
    print( paste( 'No dataset specified, defaulting to', dataset))
  }
  
  # name months 1:12 for extracting from raster
  names.months <- paste0( year.in, '-',
                          formatC( 1:12, width = 2, flag = '0'), '-',
                          '01')
  
  # extract monthly dates using function from hyspdisp
  raster.sub <- brick( subset_nc_date(  hpbl_brick = raster.in,
                                        vardate = names.months))
  
  #NARR dataset requires rotating
  if( dataset != 'NARR')
    raster.sub <- rotate( raster.sub)
  
  # take annual mean
  if( avg.period == 'year'){
    out <- stackApply( raster.sub, indices = rep( 1, 12), fun = mean)
  } else
    out <- raster.sub
  
  return( out)
}

#======================================================================#
## functions to get meteorology data
# trim data over US, create raster object
#======================================================================#
usa.functioner <- function( year.in = 2005,
                            list.met,
                            dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR'),
                            avg.period = 'year',
                            return.usa.mask = F,
                            return.usa.sub = T){
  
  # extract year
  if( avg.period == 'year'){
    mets <- lapply( list.met,
                    extract_year.fn,
                    year.in = year.in,
                    avg.period = avg.period,
                    dataset = dataset)
  } else
    mets <- lapply( list.met,
                    extract_year.fn,
                    year.in = year.in,
                    avg.period = avg.period,
                    dataset = dataset)
  
  crs.str <- projection( list.met[[1]])
  crs.usa <- crs( crs.str)
  
  # convert temp to celcius
  mets$temp <- mets$temp - 273.15
  
  # calculate windspeed
  # calculate meteorology wind angle (0 is north wind)
  # http://weatherclasses.com/uploads/3/6/2/3/36231461/computing_wind_direction_and_speed_from_u_and_v.pdf
  if( 'uwnd' %in% names( list.met) & 'vwnd' %in% names( list.met)){
    mets$wspd <- sqrt( mets$uwnd ^ 2 + mets$vwnd ^ 2)
    mets$phi <- atan2( mets$uwnd, mets$vwnd) * 180 / pi + 180
    names( mets$wspd) <- names(mets$vwnd)
    names( mets$phi) <- names(mets$vwnd)
  }
  
  # download USA polygon from rnaturalearth
  # download USA polygon from rnaturalearth
  us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
  us_states <- st_transform( USAboundaries::us_states(), crs.str)
  us_states.sp <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]
  
  if( return.usa.mask){
    return( us_states.sp)
  }
  
  if( return.usa.sub){
    mets_crop <- rasterize( us_states.sp, mets[[1]], getCover=TRUE)
    mets_crop[mets_crop==0] <- NA
    mets_crop[!is.na( mets_crop)] <- 1
    
    # crop to USA
    if( avg.period == 'year'){
      mets.out.l <- lapply( mets, function( x){
        trim( mask(x, mets_crop, maskvalue = NA),
              padding = 1)
      })
      
      mets.out <- brick( mets.out.l)
      
    } else{ 
      names.months <- names( mets[[1]])
      mets.out.l <- lapply( mets, 
                            function( x){
                              lapply( names.months, 
                                      function( y){
                                        r <- subset( x, y)
                                        trim( mask(r, mets_crop, maskvalue = NA),
                                              padding = 1)
                                      })})
      
      mets.out.b <- lapply( names.months, 
                            function( n){
                              lapply( mets.out.l, function( l){
                                b <- brick( l)
                                subset( b, n)
                              })
                            })
      
      # each month is a brick
      mets.out <- lapply( mets.out.b, brick)
      names( mets.out) <- names.months
    }
    
    return( mets.out)
  } else{
    if( avg.period == 'year'){
      mets.out <- brick( mets)
    } else{ 
      
      # reorganize - list of months
      names.months <- names( mets[[1]])
      mets.out <- lapply( names.months, 
                          function( name.ext, X){
                            brick( lapply( X, '[[', name.ext))
                          }, mets)
      names( mets.out) <- names.months
    }
    return( mets.out)
  }
}

#======================================================================#
# define the spBayes model function
#======================================================================#
splM.hyads.ddm <- function( dummy.n = 1,
                            coords = as.matrix( hyads2005.dt[,.( x, y)]),
                            Y = ddm2005.dt$X2005,
                            X = data.table( intercept = 1, hyads = hyads2005.dt$X2005),
                            seed.n = NULL,
                            ...){
  
  set.seed( seed.n)
  
  quants <- function(x){
    quantile(x, prob=c(0.5, 0.025, 0.975))
  }
  
  # holdout fraction
  ho.frac <- .1
  
  # define number of samples
  n.samples <- 5000
  
  # define priors
  starting <- list("tau.sq"=1, "sigma.sq"=1, "phi"=6)
  tuning <- list("tau.sq"=0.01, "sigma.sq"=0.01, "phi"=0.1)
  priors <- list("beta.Flat", "tau.sq.IG"=c(2, 1),
                 "sigma.sq.IG"=c(2, 1), "phi.Unif"=c(3, 30))
  
  # define holdout parameters
  ho <- sample( 1:length( Y), ceiling( ho.frac * length( Y)))
  
  #convert X & Y to matrices
  X.m <- as.matrix( X)
  Y.m <- as.matrix( Y)
  
  # define inputs
  coords.ho <- coords[ho,]
  coords.tr <- coords[-ho,]
  Y.ho <- Y.m[ho]
  Y.tr <- Y.m[-ho]
  X.ho <- X.m[ho,]
  X.tr <- X.m[-ho,]
  
  # define burn in
  burn.in <- floor(0.75*n.samples)
  
  # train the model
  m.i <- spLM( Y.tr ~ X.tr - 1, coords = coords.tr, 
               modified.pp = TRUE, ...,
               starting = starting, tuning = tuning, priors = priors, 
               cov.model = "exponential",
               n.samples = n.samples, n.report = 2500)
  
  # recover estimates of beta and theta
  rt.rec <- system.time( {
    m.i.rec <- spRecover(m.i, start=burn.in, thin=5, n.report=100)
  })
  
  # return simulated y.hat's
  rt.y.hat <- system.time( {
    m.i.pred <- spPredict( m.i, start=burn.in, thin=2, pred.covars = X.ho, 
                           pred.coords=coords.ho, verbose=FALSE)
  })
  
  # find estimates of beta, theta, w.hat, and y.hat
  beta.hat <- round(summary(m.i.rec$p.beta.recover.samples)$quantiles[c(3,1,5)],6)
  theta.hat <- round( summary( window( m.i$p.theta.samples, 
                                       start = burn.in))$quantiles[, c( 3, 1, 5)], 2)
  w.hat <- apply(m.i.rec$p.w.recover.samples, 1, median)
  y.hat <- data.table( t( apply(m.i.pred$p.y.predictive.samples, 1, quants)))
  
  # set up evaluation data.table
  Y.ho.tr.dt <- data.table( cbind( coords.ho, y.hat[,`50%`], Y.ho))
  setnames( Y.ho.tr.dt, 'V3', 'Y.hat')
  
  # rasterize output for plots
  w.hat.raster <- rasterFromXYZ( cbind( coords.tr, w.hat))
  y.hat.raster <- rasterFromXYZ( Y.ho.tr.dt[, .( x, y, Y.hat)])
  Y.tr.raster <- rasterFromXYZ( cbind( coords.tr, Y.tr))
  Y.ho.raster <- rasterFromXYZ( Y.ho.tr.dt[, .( x, y, Y.ho)])
  
  # mcmc plots
  # plot( m.i$p.theta.samples)
  # plot( m.i.rec$p.beta.recover.samples)
  
  # spatial areas of y.hat - Y.ho
  par(mfrow=c(1,3))
  plot( w.hat.raster, main = "w.hat (spatial adjustment term)")
  points( m.i$knot.coords, cex=1)
  plot( y.hat.raster - Y.ho.raster, main = "Y.hat - Y.ho")
  plot( (y.hat.raster - Y.ho.raster) / Y.ho.raster, main = "(Y.hat - Y.ho) / Y.ho")
  par(mfrow=c(1,1))
  
  # check out simpler models - set up inputs
  names.covars <- names( X)
  XY.tr <- data.table( Y = Y.tr, X.tr)
  XY.ho <- data.table( Y = Y.ho, X.ho)
  setnames( XY.tr, names(XY.tr)[names(XY.tr) %ni% 'Y'], names.covars)
  setnames( XY.ho, names(XY.tr)[names(XY.tr) %ni% 'Y'], names.covars)
  
  # check out simpler models - define them
  form <- as.formula( paste( 'Y ~ -1 +', paste( names.covars, collapse = '+')))
  m.lm <- lm( form, data = XY.tr)
  m.mean <- mean( Y.tr / XY.tr$hyads)
  
  # check out simpler models - get predicted Y.ho
  y.hat.lm <- predict( m.lm, newdata = XY.ho)
  y.hat.mean <- XY.ho$hyads * m.mean
  
  # calculate evaluation metrics
  eval.fn <- function( Yhat, Yact, mod.name){
    num.diff <- sum( Yhat - Yact)
    abs.diff <- sum( abs( Yhat - Yact))
    denom <- sum( Yact)
    metrics <- data.table( mod.name = mod.name,
                           NMB = num.diff / denom,
                           NME = abs.diff / denom,
                           MB   = num.diff / length( ho),
                           RMSE = sqrt( sum( ( Yhat - Yact) ^ 2) / length( Y.ho)),
                           R = cor( Yhat, Yact))
    return( metrics)
  }
  metrics.out <- rbind( eval.fn( Y.ho.tr.dt$Y.hat, Y.ho, 'spLM'),
                        eval.fn( y.hat.lm, Y.ho, 'lm'),
                        eval.fn( y.hat.mean, Y.ho, 'mean'))
  
  out <- list( metrics = metrics.out,
               runtimes = data.table( ncells.tr = length( Y.tr),
                                      train = round(m.i$run.time[3]/60,3),
                                      recover = round(rt.rec[3]/60,3),
                                      est.y.hat = round(rt.y.hat[3]/60,3)))
  
  return( out)
}

#======================================================================#
# define the linear model holdout function
#======================================================================#
lm.hyads.ddm.holdout <- function( seed.n = NULL,
                                  dat.stack = dats2005.s.small,
                                  dat.stack.pred = NULL,
                                  y.name = 'cmaq.ddm',
                                  x.name = 'hyads',
                                  name.idwe = 'tot.sum',
                                  covars.names = NULL, #c( 'temp', 'apcp'),
                                  ho.frac = .1,
                                  return.mods = F,
                                  ...){
  # define eval function
  eval.fn <- function( Yhat, Yact, mod.name){
    num.diff <- sum( Yhat - Yact, na.rm = T)
    abs.diff <- sum( abs( Yhat - Yact), na.rm = T)
    denom <- sum( Yact, na.rm = T)
    metrics <- data.table( mod.name = mod.name,
                           NMB = num.diff / denom,
                           NME = abs.diff / denom,
                           MB   = num.diff / length( Yhat),
                           RMSE = sqrt( sum( ( Yhat - Yact) ^ 2, na.rm = T) / length( Yhat)),
                           R = cor( Yhat, Yact, use = 'complete.obs'),
                           R.s = cor( Yhat, Yact, use = 'complete.obs', method = 'spearman'))
    return( metrics)
  }
  
  set.seed( seed.n)
  
  # if no covar names provided, use all covariates that aren't x or y
  if( is.null( covars.names))
    covars.names <- names( dat.stack)[names( dat.stack) %ni% c( x.name, y.name)]
  
  # define holdout parameters
  N <- ncell( dat.stack)
  ho <- sample( 1:N, ceiling( ho.frac * N))
  
  # extract coordinates
  dat.coords <- coordinates( dat.stack)
  dat.coords.ho <- data.table( dat.coords[ho,])
  dat.coords.tr <- data.table( dat.coords[-ho,])
  
  # define inputs
  dat.stack.ho <- data.table( values( dat.stack)[ ho,])
  dat.stack.tr <- data.table( values( dat.stack)[-ho,])
  
  # special case for ho is zero
  if( length( ho) == 0){
    # extract coordinates
    dat.coords.ho <- data.table( dat.coords)
    dat.coords.tr <- data.table( dat.coords)
    
    # define inputs
    dat.stack.ho <- data.table( dat.coords.ho, values( dat.stack))
    dat.stack.tr <- data.table( dat.coords.tr, values( dat.stack))
  } 
  if( !is.null( dat.stack.pred)){
    # extract coordinates
    dat.coords.ho <- data.table( coordinates( dat.stack.pred))
    
    # define inputs
    dat.stack.ho <- data.table( dat.coords.ho, values( dat.stack.pred))
  } 
  
  # create log variables
  dat.stack.tr[, y.name.log := log( get( y.name))]
  dat.stack.ho[, y.name.log := log( get( y.name))]
  
  # check out linear regression models - define them
  form.ncv <- as.formula( paste( 'y.name.log', '~', x.name))
  form.cv_single_poly <- 
    as.formula( paste( 'y.name.log', '~ poly(', x.name, ', 2) +', x.name, ': (', 
                       paste( c( covars.names), 
                              collapse = '+'), ')')) #, '+ s( x, y, k = 20)'
  form.cv_single <-  
    as.formula( paste( 'y.name.log', '~ ', x.name, '+', x.name, ': (', 
                       paste( c( covars.names), 
                              collapse = '+'), ')'))
  form.cv_five <- 
    as.formula( paste( 'y.name.log', '~ ', x.name, '+', x.name, ': (', 
                       paste( c( covars.names), 
                              collapse = '+'), ')^5'))
  
  # train the models
  lm.ncv <- lm( form.ncv, data = dat.stack.tr)
  lm.cv_single <-  lm( form.cv_single,  data = dat.stack.tr)
  lm.cv_single_poly <-  lm( form.cv_single_poly,  data = na.omit( dat.stack.tr))
  lm.cv_five <-  lm( form.cv_five,  data = dat.stack.tr)
  
  # check out simpler models - get predicted Y.ho
  y.ho <- unlist( dat.stack.ho[,..y.name])
  y.hat.lm.ncv  <- predict( lm.ncv,  type = 'response',  newdata = dat.stack.ho, se.fit = T)
  y.hat.lm.cv_single <- predict( lm.cv_single, type = 'response', newdata = dat.stack.ho, se.fit = T)
  y.hat.lm.cv_single_poly <- predict( lm.cv_single_poly,  type = 'response', newdata = dat.stack.ho, se.fit = T)
  y.hat.lm.cv_five <- predict( lm.cv_five,  type = 'response', newdata = dat.stack.ho, se.fit = T)
  
  # set up evaluation data.table
  Y.ho.hat <- data.table( dat.coords.ho, y.ho, 
                          y.hat.lm.ncv = y.hat.lm.ncv$fit, 
                          y.hat.lm.cv_single = y.hat.lm.cv_single$fit, 
                          y.hat.lm.cv_single_poly = y.hat.lm.cv_single_poly$fit, 
                          y.hat.lm.cv_five = y.hat.lm.cv_five$fit)
  Y.ho.hat.bias <- data.table( dat.coords.ho, y.ho, 
                               y.hat.lm.ncv = y.hat.lm.ncv$fit - y.ho, 
                               y.hat.lm.cv_single = y.hat.lm.cv_single$fit - y.ho, 
                               y.hat.lm.cv_single_poly = y.hat.lm.cv_single_poly$fit - y.ho, 
                               y.hat.lm.cv_five = y.hat.lm.cv_five$fit - y.ho)
  Y.ho.hat.se <- data.table( dat.coords.ho, 
                             y.hat.lm.ncv = y.hat.lm.ncv$se.fit, 
                             y.hat.lm.cv_single = y.hat.lm.cv_single$se.fit,  
                             y.hat.lm.cv_single_poly = y.hat.lm.cv_single_poly$se.fit,  
                             y.hat.lm.cv_five = y.hat.lm.cv_five$se.fit)
  
  # rasterize output for plots
  crs.in <- crs( dat.stack)
  Y.ho.hat.raster <- projectRaster( rasterFromXYZ( Y.ho.hat, crs = crs.in), dat.stack)
  Y.ho.hat.se.raster <- projectRaster( rasterFromXYZ( Y.ho.hat.se, crs = crs.in), dat.stack)
  Y.ho.hat.bias.raster <- projectRaster( rasterFromXYZ( Y.ho.hat.bias, crs = crs.in), dat.stack)
  
  # calculate evaluation metrics
  metrics.out <- rbind( eval.fn( exp( Y.ho.hat$y.hat.lm.ncv), y.ho, 'lm.ncv'),
                        eval.fn( exp( Y.ho.hat$y.hat.lm.cv_single), y.ho, 'lm.cv_single'),
                        eval.fn( exp( Y.ho.hat$y.hat.lm.cv_single_poly), y.ho, 'lm.cv_single_poly'),
                        eval.fn( exp( Y.ho.hat$y.hat.lm.cv_five), y.ho, 'lm.cv_five'))
  
  # listify the models
  if( return.mods)
    out <- list( metrics = metrics.out, 
                 model.lm.ncv = lm.ncv,
                 model.lm.cv_single = lm.cv_single,
                 model.lm.cv_single_poly = lm.cv_single_poly,
                 model.lm.cv_five = lm.cv_five,
                 Y.ho.hat.raster = Y.ho.hat.raster, 
                 Y.ho.hat.se.raster = Y.ho.hat.se.raster,
                 Y.ho.hat.bias.raster = Y.ho.hat.bias.raster)
  if( !return.mods)
    out <- list( metrics = metrics.out, 
                 Y.ho.hat.raster = Y.ho.hat.raster, 
                 Y.ho.hat.se.raster = Y.ho.hat.se.raster,
                 Y.ho.hat.bias.raster = Y.ho.hat.bias.raster)
  
  return( out)
}

#======================================================================#
# extract a mean from a list of lists of rasters
#======================================================================#
mean.lol <- function( lol, layer1, layer2, plot.out = TRUE){
  first.lol <- lapply( lol, '[[', layer1)
  second.lol <- lapply( first.lol, '[[', layer2)
  mean.lol <- mean( stack( second.lol), na.rm = T)
  
  if( plot.out)
    plot( mean.lol, main = layer2)
  return( mean.lol)
}

#======================================================================#
# ggplot a raster
#======================================================================#
ggplot.a.raster <- function( ..., bounds = NULL, facet.names = NULL, mask.raster = NULL, 
                             nrow. = NULL, ncol. = NULL, legend.name = NULL, theme.obj = theme()){
  
  in.x <- list( ...)
  
  if( length( in.x) == 1)
    in.x <- in.x[[1]]
  
  if( is.null( facet.names))
    facet.names <- lapply( in.x, names)
  
  if( is.null( names( in.x)) & !is.null( facet.names))
    names( in.x) <- facet.names
  
  in.x.crop <- in.x
  if( !is.null( mask.raster) & is.list( in.x))
    in.x.crop <- lapply( in.x, function( X, mask.raster.){
      X.mask <- mask( X, mask.raster.)
      X.crop <- crop( X.mask, mask.raster.)
      return( X.crop)
    }, mask.raster)
  if( !is.null( mask.raster) & !is.list( in.x)){
    in.x.mask <- mask( in.x, mask.raster)
    in.x.crop <- crop( in.x.mask, mask.raster)
  }
  
  dat.dt <- rbindlist( lapply( facet.names, function( x.name, x.list) {
    x <- x.list[[x.name]]
    r_points <- rasterToPoints( x)
    r_dt <- data.table( r_points)[, name.in := x.name]
    setnames( r_dt, names( x), 'z')
    return( r_dt)
  }, in.x.crop))
  
  if( !is.null( facet.names))
    dat.dt[, name.in := factor( name.in, levels = facet.names)]
  
  ggplot( dat.dt) +
    geom_tile( aes( x = x, y = y, fill = z)) + 
    scale_fill_viridis( name = legend.name, limits = bounds, oob = scales::squish) + 
    facet_wrap( . ~ name.in, nrow = nrow., ncol = ncol.) + 
    # expand_limits( fill = 0) +
    theme_bw() + 
    theme( axis.text = element_blank(),
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           legend.position = 'bottom',
           legend.key.width = unit( ncol( in.x[[1]])/3000, 'npc'),
           panel.grid = element_blank(),
           strip.background = element_blank()) + 
    theme.obj
}

#======================================================================#
# do the predictions for many months
#======================================================================#
month.trainer <- function( name.m = names( mets2005.m)[1], 
                           name.p = names( mets2006.m)[1], 
                           name.x,
                           y.m, 
                           ddm.m = ddm.m.all, 
                           mets.m = mets.m.all,
                           emiss.m = d_nonegu.r,
                           idwe.m = idwe.m,
                           .mask.use = NULL,
                           cov.names = c( "temp", "rhum", "vwnd", "uwnd", "wspd")){ #, names( d_nonegu.r))){
  
  # create training dataset
  ddm.use <- ddm.m[[name.m]]
  hyads.use <- y.m[[name.m]]
  mets.use <- mets.m[[name.m]]
  emiss.use <- emiss.m
  idwe.use <- idwe.m[[name.m]]
  
  # create prediction dataset
  ddm.use.p <- ddm.m[[name.p]]
  hyads.use.p <- y.m[[name.p]]
  mets.use.p <- mets.m[[name.p]]
  idwe.use.p <- idwe.m[[name.p]]
  
  # fix names
  names( ddm.use)     <- 'cmaq.ddm'
  names( ddm.use.p)   <- 'cmaq.ddm'
  names( hyads.use)   <- name.x
  names( hyads.use.p) <- name.x
  names( idwe.use)   <- 'tot.sum.idwe'
  names( idwe.use.p) <- 'tot.sum.idwe'
  
  # combine each dataset as stacks
  dat.s <- project_and_stack( hyads.use,   ddm.use,   idwe.use,   
                              mets.use, emiss.use, mask.use = .mask.use)
  dat.p <- project_and_stack( hyads.use.p, ddm.use.p, idwe.use.p, 
                              mets.use.p, emiss.use, mask.use = .mask.use)
  
  # do the modeling
  pred <- lm.hyads.ddm.holdout( dat.stack = dat.s, dat.stack.pred = dat.p, x.name = name.x,
                                ho.frac = 0, covars.names = cov.names, 
                                name.idwe = 'tot.sum.idwe', return.mods = T)
  
  return( pred)
}


#======================================================================#
# project and stack in one command
# ## need to update masking in all functions - cells with centroid not covered are cropped
##   https://gis.stackexchange.com/questions/255025/r-raster-masking-a-raster-by-polygon-also-remove-cells-partially-covered
## should probably update usa mask too - need to use USAboundaries for consistency
#======================================================================#
project_and_stack <- function( ..., mask.use = NULL){
  list.r <- list( ...)
  
  if( length( list.r) > 1)
    for( r in 2: length( list.r)){
      list.r[[r]] <- projectRaster( list.r[[r]], list.r[[1]], alignOnly = F)
    }
  
  # mask over usa
  if( !is.null( mask.use)){
    mask.use <- spTransform( mask.use, crs( list.r[[1]]))
    mask_crop <- rasterize( mask.use, list.r[[1]][[1]], getCover=TRUE)
    mask_crop[mask_crop==0] <- NA
    mask_crop[!is.na( mask_crop)] <- 1
    
    list.out <- lapply( list.r, function( x){
      x1 <- reclassify( x, c( NA, NA, 0))
      # x[is.na(x)] <- 0
      trim( mask(x1, mask_crop, maskvalue = NA),
            padding = 1)
    })
    
  } else
    list.out <- list.r
  
  return( stack( list.out))
}

#======================================================================#
## function for converting individual unit concentrations to ugm3
# inputs - filename of grid w/ unit or summed impacts
#        - model
#        - covariate raster
#        - 'x' name in the raster from model
# 
# output - data table of state, pop-wgted impact for all input units
#======================================================================#
state_exposurer <- function( 
    month.n,
    fstart,
    year.m = 2006,
    model.dataset = preds.mon.idwe06w05,
    model.name = 'model.cv', #'model.gam'
    name.x = 'idwe',
    mask.use = mask.usa,
    ddm.m = ddm.m.all,
    mets.m = mets.m.all,
    emiss.m = d_nonegu.r,
    idwe.m. = idwe.m,
    hyads.m = hyads.m.all,
    grid_pop.r = grid_popwgt.r,
    state_pops = copy( us_states.pop.dt),
    p4s,
    take.diff = F,
    xboost = F
){
  message( paste( 'Converting', month.name[month.n]))
  
  month.N <- formatC( month.n, width = 2, flag = '0')
  name.m <- paste0( 'X', year.m, '.', month.N, '.01')
  model.m <- paste0( 'X2005.', month.N, '.01')
  popyr.name <- paste0( 'X', year.m)
  name.Date <- as.Date( name.m, format = 'X%Y.%m.%d')
  
  if( name.x == 'idwe'){
    fname <- paste0( fstart, year.m, '_', month.n, '.csv')
    name.dat <- 'idwe'
  } else{
    fname <- paste0( fstart, year.m, '_', month.N, '.csv')
    name.dat <- 'hyads'
  }
  
  # create prediction dataset
  # ddm.use.p <- ddm.m[[name.m]]
  mets.use.p <- mets.m[[name.m]]
  idwe.use.p <- idwe.m.[[name.m]]
  hyads.use.p <- hyads.m[[name.m]]
  
  # fix names
  # names( ddm.use.p)   <- 'cmaq.ddm'
  names( hyads.use.p) <- 'hyads'
  names( idwe.use.p) <- 'idwe'
  
  # rename population raster
  grid_pop.r <- grid_pop.r[[popyr.name]]
  names( grid_pop.r) <- 'pop'
  
  dat.s <- project_and_stack( #ddm.use.p,   
    hyads.use.p,   idwe.use.p,   
    mets.use.p, emiss.m, grid_pop.r, mask.use = mask.use)
  
  # assign state names to raster
  mask.r <- rasterize( mask.use[,'state_abbr'], dat.s[[1]])
  mask.a <- levels( mask.r)[[1]]
  dat.s$ID <- mask.r
  
  # read in, remove 
  x.in1 <- fread( fname, drop = c( 'V1'))
  x.in1[is.na( x.in1)] <- 0
  suppressWarnings( x.in1[, `:=` ( yearmon = NULL, yearmonth = NULL)])
  x.in2 <- x.in1[, colSums(x.in1) != 0, with = F]
  suppressWarnings( x.in2[, `:=` ( x = NULL, y = NULL)])
  x.in <- cbind( x.in1[, .( x, y)], x.in2)
  
  # rasterize new x file
  x.r <- rasterFromXYZ( x.in, crs = p4s)
  x.n <- paste0( 'X', names( x.in)[!(names( x.in) %in% c( 'x', 'y'))])
  x.n <- gsub( '^XX', 'X', x.n)
  names( x.r) <- x.n
  x.proj <-  project_and_stack( dat.s[[1]], x.r, mask.use = mask.use)
  names( x.proj)[2:dim( x.proj)[3]] <- x.n
  
  # pick out the appropriate model
  model.use <- model.dataset[ model.name, model.m][[1]]
  
  # predict the base scenario (zero hyads/idwe)
  dat.use0<- copy( dat.s)
  dat.use0[[name.dat]] <- 0
  dat.coords <- coordinates( dat.s)
  dat_raw0.dt <- data.table( cbind( dat.coords, values( dat.use0)))
  
  # xboost requires special treatment
  if( xboost){
    dat_raw0.dt.trim <- dat_raw0.dt[, model.use$feature_names, with = F]
    xhold1c <- xgb.DMatrix( as.matrix( dat_raw0.dt.trim))
    dat.pred0 <- predict( model.use, newdata = xhold1c)
  } else
    dat.pred0 <- predict( model.use, newdata = dat_raw0.dt)
  
  dats0.r <- rasterFromXYZ( data.table( dat.coords, dat.pred0), crs = p4s)
  
  # do the predictions
  pred_pm.r <- brick( pbmcapply::pbmclapply( x.n, function( n){ 
    gc()
    # assign unit to prediction dataset
    dat.use <- copy( dat.s)
    dat.use[[name.dat]] <- x.proj[[n]]
    
    # if zero impacts, return raster with only zeros
    if( sum( values(x.proj[[n]]), na.rm = T) == 0)
      return( x.proj[[n]])
    
    # set up the dataset
    dat_raw.dt <- data.table( cbind( dat.coords, values( dat.use)))
    
    # do the predictions
    if( xboost){
      dat_raw.dt.trim <- dat_raw.dt[, model.use$feature_names, with = F]
      xhold.pred <- xgb.DMatrix( as.matrix( dat_raw.dt.trim))
      dat.pred <- predict( model.use, newdata = xhold.pred)
    } else
      dat.pred <- predict( model.use, newdata = dat_raw.dt)
    
    # rasterize
    dats.r <- rasterFromXYZ( data.table( dat.coords, dat.pred), crs = p4s)
    
    #take difference from base
    if( take.diff){
      dats.r2 <- dats.r - dats0.r
    } else 
      dats.r2 <- dats.r
    
    names( dats.r2) <- n
    return( dats.r2)
  }))
  
  # multiply by population
  pred_popwgt.r <- pred_pm.r * dat.s[['pop']]
  names( pred_popwgt.r) <- names( pred_pm.r)
  
  #calculate raw average for the entire domain
  pred_pm.us <- colMeans( data.table( values( pred_pm.r)), na.rm = T)
  pred_pm.us.dt <- data.table( uID = names( pred_pm.us),
                               mean_pm = pred_pm.us,
                               state_abbr = 'US',
                               ID = 100)
  
  #calculate pop-weightedverage for the entire domain
  pred_pm.pw.us <- colSums( na.omit( data.table( values( pred_popwgt.r)), na.rm = T))
  pred_pm.pw.us.dt <- data.table( uID = names( pred_pm.pw.us),
                                  mean_popwgt = pred_pm.pw.us,
                                  state_abbr = 'US',
                                  ID = 100)
  
  # calculate raw average by state
  pred_pm.r$ID <- mask.r
  pred_pm.dt <- data.table( values( pred_pm.r))
  pred_pm.dt.m <- na.omit( melt( pred_pm.dt, id.vars = 'ID', variable.name = 'uID'))
  pred_pm.dt.s <- pred_pm.dt.m[, .( mean_pm = mean( value)), by = .( ID, uID)]
  pred_pm.dt.s <- merge( pred_pm.dt.s, mask.a, by = 'ID')
  
  # calculate population-weighted by state
  pred_popwgt.r$ID <- mask.r
  pred_popwgt.dt <- data.table( values( pred_popwgt.r))
  pred_popwgt.dt.m <- na.omit( melt( pred_popwgt.dt, id.vars = 'ID', variable.name = 'uID'))
  pred_popwgt.dt.s <- pred_popwgt.dt.m[, .( mean_popwgt = sum( value)), by = .( ID, uID)]
  pred_popwgt.dt.s <- merge( pred_popwgt.dt.s, mask.a, by = 'ID')
  
  # bind with united states pops
  pred_pm <- rbind( pred_pm.dt.s, pred_pm.us.dt)
  pred_pm.pw <- rbind( pred_popwgt.dt.s, pred_pm.pw.us.dt)
  
  # now just divide by each state's total population
  setnames( state_pops, popyr.name, 'pop_amnt')
  state_pops_lite <- state_pops[, .( state_abbr, pop_amnt)]
  pop.tot <- data.table( state_abbr = 'US', pop_amnt = sum( state_pops$pop_amnt))
  state_pops_lite <- rbind( state_pops_lite, pop.tot)
  pred_popwgt.out <- merge( pred_pm.pw, state_pops_lite, by = 'state_abbr')
  
  # divide by total population
  pred_popwgt.out[, `:=` ( popwgt = mean_popwgt / pop_amnt,
                           month = name.Date)]
  
  # merge pop-weighted and raw dataset
  out <- merge( pred_popwgt.out, pred_pm, by = c( 'state_abbr', 'ID', 'uID'))
  
  return( list( popwgt_states = out, pred_pm.r = pred_pm.r, zero_out.r = dats0.r))
  
}

#======================================================================#
## same as above, but for a year
#======================================================================#
state_exposurer.year <- function( 
    fname,
    year.m = 2006,
    model.use = preds.ann.hyads06w05$model.gam,
    name.x = 'idwe',
    mask.use = mask.usa,
    dat.a = dats2006.a,
    grid_pop.r = grid_popwgt.r,
    state_pops = copy( us_states.pop.dt),
    p4s,
    take.diff = F,
    xboost = F,
    raw = F
){
  message( paste( 'Converting', year.m))
  
  # rename population raster
  popyr.name <- paste0( 'X', year.m)
  grid_pop.r <- grid_pop.r[[popyr.name]]
  names( grid_pop.r) <- 'pop'
  
  # assign state names to raster
  mask.r <- rasterize( mask.use[,'state_abbr'], dat.a[[1]])
  mask.a <- levels( mask.r)[[1]]
  dat.a$ID <- mask.r
  dat.a <- project_and_stack( dat.a, grid_pop.r)
  
  # read in, rasterize new x file
  if( name.x == 'hyads'){
    x.in <- fread( fname, drop = c( 'V1', 'year.E', 'year.H'))
    x.cast <- dcast( x.in, x + y ~ uID, value.var = 'hyads')
    name.dat <- 'hyads'
  }
  if( name.x == 'idwe'){
    x.cast <- fread( fname, drop = c( 'V1'))
    name.dat <- 'idwe'
  }
  
  # cast and rasterize
  x.r <- rasterFromXYZ( x.cast, crs = p4s)
  x.n <- paste0( 'X', names( x.cast)[!(names( x.cast) %in% c( 'x', 'y'))])
  x.n <- gsub( '^XX', 'X', x.n)
  names( x.r) <- x.n
  x.proj <-  project_and_stack( dat.a[[1]], x.r, mask.use = mask.use)
  
  # predict the base scenario (zero hyads/idwe)
  dat.use0<- copy( dat.a)
  dat.use0[[name.dat]] <- 0
  dat.coords <- coordinates( dat.a)
  dat_raw0.dt <- data.table( cbind( dat.coords, values( dat.use0)))
  
  # do the prediction if not taking raw values
  if( !raw){
    # xboost requires special treatment
    if( xboost){
      dat_raw0.dt.trim <- dat_raw0.dt[, model.use$feature_names, with = F]
      xhold1c <- xgb.DMatrix( as.matrix( dat_raw0.dt.trim))
      dat.pred0 <- predict( model.use, newdata = xhold1c)
    } else
      dat.pred0 <- predict( model.use, newdata = dat_raw0.dt)
    
    dats0.r <- rasterFromXYZ( data.table( dat.coords, dat.pred0), crs = p4s)
  } else{
    dats0.r <- rasterFromXYZ( dat_raw0.dt[, c( 'x', 'y', name.dat), with = F], crs = p4s)
  }
  
  
  # do the predictions
  pred_pm.r <- brick( pbmcapply::pbmclapply( x.n, function( n){ 
    gc()
    # if zero impacts, return raster with only zeros
    if( sum( values(x.proj[[n]]), na.rm = T) == 0 | raw)
      return( x.proj[[n]])
    
    # assign unit to prediction dataset
    dat.use <- copy( dat.a)
    dat.use[[name.dat]] <- x.proj[[n]]
    
    # set up the dataset
    dat_raw.dt <- data.table( cbind( dat.coords, values( dat.use)))
    
    # do the predictions
    if( xboost){
      dat_raw.dt.trim <- dat_raw.dt[, model.use$feature_names, with = F]
      xhold.pred <- xgb.DMatrix( as.matrix( dat_raw.dt.trim))
      dat.pred <- predict( model.use, newdata = xhold.pred)
    } else
      dat.pred <- predict( model.use, newdata = dat_raw.dt)
    
    # rasterize
    dats.r <- rasterFromXYZ( data.table( dat.coords, dat.pred), crs = p4s)
    
    #take difference from base
    if( take.diff){
      dats.r2 <- dats.r - dats0.r
    } else 
      dats.r2 <- dats.r
    
    names( dats.r2) <- n
    return( dats.r2)
  }))
  
  # multiply by population
  pred_popwgt.r <- pred_pm.r * dat.a[['pop']]
  names( pred_popwgt.r) <- names( pred_pm.r)
  
  #calculate raw average for the entire domain
  pred_pm.us <- colMeans( data.table( values( pred_pm.r)), na.rm = T)
  pred_pm.us.dt <- data.table( uID = names( pred_pm.us),
                               mean_pm = pred_pm.us,
                               state_abbr = 'US',
                               ID = 100)
  
  #calculate pop-weightedverage for the entire domain
  pred_pm.pw.us <- colSums( na.omit( data.table( values( pred_popwgt.r)), na.rm = T))
  pred_pm.pw.us.dt <- data.table( uID = names( pred_pm.pw.us),
                                  mean_popwgt = pred_pm.pw.us,
                                  state_abbr = 'US',
                                  ID = 100)
  
  # calculate raw average by state
  pred_pm.r$ID <- mask.r
  pred_pm.dt <- data.table( values( pred_pm.r))
  pred_pm.dt.m <- na.omit( melt( pred_pm.dt, id.vars = 'ID', variable.name = 'uID'))
  pred_pm.dt.s <- pred_pm.dt.m[, .( mean_pm = mean( value)), by = .( ID, uID)]
  pred_pm.dt.s <- merge( pred_pm.dt.s, mask.a, by = 'ID')
  
  # calculate population-weighted by state
  pred_popwgt.r$ID <- mask.r
  pred_popwgt.dt <- data.table( values( pred_popwgt.r))
  pred_popwgt.dt.m <- na.omit( melt( pred_popwgt.dt, id.vars = 'ID', variable.name = 'uID'))
  pred_popwgt.dt.s <- pred_popwgt.dt.m[, .( mean_popwgt = sum( value)), by = .( ID, uID)]
  pred_popwgt.dt.s <- merge( pred_popwgt.dt.s, mask.a, by = 'ID')
  
  # bind with united states pops
  pred_pm <- rbind( pred_pm.dt.s, pred_pm.us.dt)
  pred_pm.pw <- rbind( pred_popwgt.dt.s, pred_pm.pw.us.dt)
  
  # now just divide by each state's total population
  setnames( state_pops, popyr.name, 'pop_amnt')
  state_pops_lite <- state_pops[, .( state_abbr, pop_amnt)]
  pop.tot <- data.table( state_abbr = 'US', pop_amnt = sum( state_pops$pop_amnt))
  state_pops_lite <- rbind( state_pops_lite, pop.tot)
  pred_popwgt.out <- merge( pred_pm.pw, state_pops_lite, by = 'state_abbr')
  
  # divide by total population
  pred_popwgt.out[, `:=` (popwgt = mean_popwgt / pop_amnt,
                          year = year.m)]
  
  # merge pop-weighted and raw dataset
  out <- merge( pred_popwgt.out, pred_pm, by = c( 'state_abbr', 'ID', 'uID'))
  
  return( list( popwgt_states = out, pred_pm.r = pred_pm.r, zero_out.r = dats0.r))
}

#======================================================================#
## calculate evaluation metrics
#======================================================================#
evals.fn <- function( Yhat, Yact){
  num.diff <- sum( Yhat - Yact)
  abs.diff <- sum( abs( Yhat - Yact))
  denom <- sum( Yact)
  metrics <- data.table( NMB = num.diff / denom,
                         NME = abs.diff / denom,
                         MB   = num.diff / length( Yhat),
                         ME   = abs.diff / length( Yhat),
                         RMSE = sqrt( sum( ( Yhat - Yact) ^ 2) / length( Yhat)),
                         R.p = cor( Yhat, Yact, method = 'pearson'),
                         R.s = cor( Yhat, Yact, method = 'spearman'))
  return( metrics)
}

#======================================================================#
## function for converting individual unit concentrations to ugm3
# inputs - filename of grid w/ unit or summed impacts
#        - model
#        - covariate raster
#        - 'x' name in the raster from model
# 
# output - data table of state, pop-wgted impact for all input units
#======================================================================#
hyads_to_pm25 <- function( 
    month.n = NULL,
    fstart,
    year.m = 2006,
    model.dataset = preds.mon.idwe06w05,
    model.name = "model.cv", 
    name.x = 'idwe',
    mask.use = mask.usa,
    ddm.m = ddm.m.all,
    mets.m = mets.m.all,
    emiss.m = d_nonegu.r,
    idwe.m. = idwe.m,
    hyads.m = hyads.m.all,
    take.diff = T,
    p4s
){
  message( paste( 'Converting', month.name[month.n]))
  
  month.N <- formatC( month.n, width = 2, flag = '0')
  name.m <- paste0( 'X', year.m, '.', month.N, '.01')
  model.m <- paste0( 'X2005.', month.N, '.01')
  popyr.name <- paste0( 'X', year.m)
  name.Date <- as.Date( name.m, format = 'X%Y.%m.%d')
  
  if( name.x == 'idwe'){
    fname <- paste0( fstart, year.m, '_', month.n, '.csv')
    name.dat <- 'idwe'
  } else{
    fname <- paste0( fstart, year.m, '_', month.N, '.csv')
    name.dat <- 'hyads'
  }
  
  # create prediction dataset
  # ddm.use.p <- ddm.m[[name.m]]
  mets.use.p <- mets.m[[name.m]]
  idwe.use.p <- idwe.m.[[name.m]]
  hyads.use.p <- hyads.m[[name.m]]
  
  # fix names
  # names( ddm.use.p)   <- 'cmaq.ddm'
  names( hyads.use.p) <- 'hyads'
  names( idwe.use.p) <- 'idwe'
  
  # rename population raster
  grid_pop.r <- grid_pop.r[[popyr.name]]
  names( grid_pop.r) <- 'pop'
  
  dat.s <- project_and_stack( #ddm.use.p,   
    hyads.use.p,   idwe.use.p,   
    mets.use.p, emiss.m, grid_pop.r, mask.use = mask.use)
  
  # assign state names to raster
  mask.r <- rasterize( mask.use[,'state_abbr'], dat.s[[1]])
  mask.a <- levels( mask.r)[[1]]
  dat.s$ID <- mask.r
  
  # read in, remove 
  x.in1 <- fread( fname, drop = c( 'V1'))
  x.in1[is.na( x.in1)] <- 0
  suppressWarnings( x.in1[, `:=` ( yearmon = NULL, yearmonth = NULL)])
  x.in2 <- x.in1[, colSums(x.in1) != 0, with = F]
  suppressWarnings( x.in2[, `:=` ( x = NULL, y = NULL)])
  x.in <- cbind( x.in1[, .( x, y)], x.in2)
  
  # rasterize new x file
  x.r <- rasterFromXYZ( x.in, crs = p4s)
  x.n <- paste0( 'X', names( x.in)[!(names( x.in) %in% c( 'x', 'y'))])
  x.n <- gsub( '^XX', 'X', x.n)
  names( x.r) <- x.n
  x.proj <-  project_and_stack( dat.s[[1]], x.r, mask.use = mask.use)
  names( x.proj)[2:dim( x.proj)[3]] <- x.n
  
  # pick out the appropriate model
  model.use <- model.dataset[ model.name, model.m][[1]]
  
  # predict the base scenario (zero hyads/idwe)
  dat.use0<- copy( dat.s)
  dat.use0[[name.dat]] <- 0
  dat.coords <- coordinates( dat.s)
  dat_raw0.dt <- data.table( cbind( dat.coords, values( dat.use0)))
  
  # xboost requires special treatment
  if( xboost){
    dat_raw0.dt.trim <- dat_raw0.dt[, model.use$feature_names, with = F]
    xhold1c <- xgb.DMatrix( as.matrix( dat_raw0.dt.trim))
    dat.pred0 <- predict( model.use, newdata = xhold1c)
  } else
    dat.pred0 <- predict( model.use, newdata = dat_raw0.dt)
  
  dats0.r <- rasterFromXYZ( data.table( dat.coords, dat.pred0), crs = p4s)
  
  # do the predictions
  pred_pm.r <- brick( pbmcapply::pbmclapply( x.n, function( n){ 
    gc()
    # assign unit to prediction dataset
    dat.use <- copy( dat.s)
    dat.use[[name.dat]] <- x.proj[[n]]
    
    # if zero impacts, return raster with only zeros
    if( sum( values(x.proj[[n]]), na.rm = T) == 0)
      return( x.proj[[n]])
    
    # set up the dataset
    dat_raw.dt <- data.table( cbind( dat.coords, values( dat.use)))
    
    # do the predictions
    if( xboost){
      xhold.pred <- xgb.DMatrix( as.matrix( dat_raw.dt.trim))
      dat.pred <- predict( model.use, newdata = xhold.pred)
    } else
      dat.pred <- predict( model.use, newdata = dat_raw.dt)
    
    # rasterize
    dats.r <- rasterFromXYZ( data.table( dat.coords, dat.pred), crs = p4s)
    
    #take difference from base
    if( take.diff){
      dats.r2 <- dats.r - dats0.r
    } else 
      dats.r2 <- dats.r
    
    names( dats.r2) <- n
    return( dats.r2)
  }))
  
  # multiply by population
  pred_popwgt.r <- pred_pm.r * dat.s[['pop']]
  names( pred_popwgt.r) <- names( pred_pm.r)
  
  #calculate raw average for the entire domain
  pred_pm.us <- colMeans( data.table( values( pred_pm.r)), na.rm = T)
  pred_pm.us.dt <- data.table( uID = names( pred_pm.us),
                               mean_pm = pred_pm.us,
                               state_abbr = 'US',
                               ID = 100)
  
  #calculate pop-weightedverage for the entire domain
  pred_pm.pw.us <- colSums( na.omit( data.table( values( pred_popwgt.r)), na.rm = T))
  pred_pm.pw.us.dt <- data.table( uID = names( pred_pm.pw.us),
                                  mean_popwgt = pred_pm.pw.us,
                                  state_abbr = 'US',
                                  ID = 100)
  
  # calculate raw average by state
  pred_pm.r$ID <- mask.r
  pred_pm.dt <- data.table( values( pred_pm.r))
  pred_pm.dt.m <- na.omit( melt( pred_pm.dt, id.vars = 'ID', variable.name = 'uID'))
  pred_pm.dt.s <- pred_pm.dt.m[, .( mean_pm = mean( value)), by = .( ID, uID)]
  pred_pm.dt.s <- merge( pred_pm.dt.s, mask.a, by = 'ID')
  
  # calculate population-weighted by state
  pred_popwgt.r$ID <- mask.r
  pred_popwgt.dt <- data.table( values( pred_popwgt.r))
  pred_popwgt.dt.m <- na.omit( melt( pred_popwgt.dt, id.vars = 'ID', variable.name = 'uID'))
  pred_popwgt.dt.s <- pred_popwgt.dt.m[, .( mean_popwgt = sum( value)), by = .( ID, uID)]
  pred_popwgt.dt.s <- merge( pred_popwgt.dt.s, mask.a, by = 'ID')
  
  # bind with united states pops
  pred_pm <- rbind( pred_pm.dt.s, pred_pm.us.dt)
  pred_pm.pw <- rbind( pred_popwgt.dt.s, pred_pm.pw.us.dt)
  
  # now just divide by each state's total population
  setnames( state_pops, popyr.name, 'pop_amnt')
  state_pops_lite <- state_pops[, .( state_abbr, pop_amnt)]
  pop.tot <- data.table( state_abbr = 'US', pop_amnt = sum( state_pops$pop_amnt))
  state_pops_lite <- rbind( state_pops_lite, pop.tot)
  pred_popwgt.out <- merge( pred_pm.pw, state_pops_lite, by = 'state_abbr')
  
  # divide by total population
  pred_popwgt.out[, `:=` ( popwgt = mean_popwgt / pop_amnt,
                           month = name.Date)]
  
  # merge pop-weighted and raw dataset
  out <- merge( pred_popwgt.out, pred_pm, by = c( 'state_abbr', 'ID', 'uID'))
  
  return( list( popwgt_states = out, pred_pm.r = pred_pm.r, zero_out.r = dats0.r))
  
}

#======================================================================#
## function for converting individual unit concentrations to ugm3
# inputs - filename of grid w/ unit or summed impacts
#        - model
#        - covariate raster
#        - 'x' name in the raster from model
# 
# output - data table of ugm3 impacts for all input units
#======================================================================#
hyads_to_pm25_unit <- function( 
    year.m = 2006,
    month.n = NULL,
    fstart = NULL,
    fstart.total,
    fstart_out,
    model.dataset = preds.mon.idwe06w05,
    model.name = 'model.cv', # 'model.gam'
    name.x = 'hyads',
    mask.use = mask.usa,
    met.dest = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met',
    total = F,
    p4s = NULL
){
  message( paste( 'Converting', month.name[month.n], year.m))
  
  #define the met layer names, do the actual downloading
  Sys.setenv(TZ='UTC')
  layer.names <- c( "air.2m.mon.mean.nc",
                    "apcp.mon.mean.nc",
                    "rhum.2m.mon.mean.nc",
                    "vwnd.10m.mon.mean.nc",
                    "uwnd.10m.mon.mean.nc")
  names( layer.names) <- c( "temp", "apcp", "rhum", "vwnd", "uwnd")
  
  # do the data downloading
  # set destination parameter to where you want the data downloaded,
  # for example, destination = '~/Desktop'
  list.met <- lapply( layer.names,
                      downloader.fn, 
                      destination = met.dest,
                      dataset = 'NARR')
  
  # annual or month-specific actions
  if( is.null( month.n)){
    # define file nanmes
    fname <- paste0( fstart, year.m, '.fst')
    fname.total <- paste0( fstart.total, year.m, '.fst')
    fname_out <- paste0( fstart_out, year.m, '.fst')
    
    # download met data
    mets.use.p <- suppressWarnings( 
      usa.functioner( year.m, list.met, dataset = 'NARR', 
                      avg.period = 'year', return.usa.sub = F)
    )
    
    # pick out the appropriate model
    model.use <- model.dataset[[ model.name]]
    
    # get the prediction crs
    model.rast <- model.dataset$Y.ho.hat.raster
    model.csr <- crs( model.rast)
    
  } else {
    # define date/month names
    month.N <- formatC( month.n, width = 2, flag = '0')
    name.m <- paste0( 'X', year.m, '.', month.N, '.01')
    model.m <- paste0( 'X2005.', month.N, '.01')
    name.Date <- as.Date( name.m, format = 'X%Y.%m.%d')
    
    # define file nanmes
    fname <- paste0( fstart, year.m, '_', month.N, '.fst')
    fname.total <- paste0( fstart.total, year.m, '_', month.N, '.fst')
    fname_out <- paste0( fstart_out, year.m, '_', month.N, '.fst')
    
    # download met data
    mets.m <- suppressWarnings( 
      usa.functioner( year.m, list.met, dataset = 'NARR', 
                      avg.period = 'month', return.usa.sub = F)
    )
    mets.use.p <- mets.m[[name.m]]
    
    # pick out the appropriate model
    model.use <- model.dataset[ model.name, model.m][[1]]
    
    # get the prediction crs
    model.rast <- model.dataset['Y.ho.hat.raster',][[model.m]]
    model.csr <- crs( model.rast)
    
  }
  
  # create prediction dataset
  mets.use.p <- projectRaster( mets.use.p, crs = model.csr)
  dat.s <- project_and_stack( model.rast, mets.use.p, mask.use = mask.use)
  
  # read in total hyads
  hyads.total.dt <- read.fst( fname.total, columns = c( 'x', 'y', 'hyads'), as.data.table = T)
  hyads.total.r <- rasterFromXYZ( hyads.total.dt, crs = p4s)
  
  # predict with total hyads
  dat.use.s<- copy( dat.s)
  hyads.total.proj <- project_and_stack( dat.use.s, hyads.total.r, mask.use = mask.use)
  
  # predict the base scenario
  dat.coords <- coordinates( hyads.total.proj)
  dat_total.dt <- data.table( cbind( dat.coords, values( hyads.total.proj)))
  dat.total.pred <- predict( model.use, newdata = dat_total.dt) %>%
    exp()
  
  # predict a 0-hyads scenario
  dat_0.dt <- dat_total.dt[, hyads := 0]
  dat.0.pred <- predict( model.use, newdata = dat_0.dt, type = 'response') %>%
    exp()
  
  # predict pm into a raster
  pred_pm.r <- rasterFromXYZ( data.table( dat.coords, dat.total.pred), crs = p4s) %>%
    projectRaster( dat.use.s)
  pred_0.r <- rasterFromXYZ( data.table( dat.coords, dat.0.pred), crs = p4s) %>%
    projectRaster( dat.use.s)
  
  # remove mean
  pred_nomean.r <- pred_pm.r - pred_0.r
  
  #read in, project hyads
  if( total){
    
    coords.out <- coordinates( pred_nomean.r)
    vals.out <- values( pred_nomean.r)
    
    # write out the data.table as fst
    pred_pm.dt <- data.table( cbind( coords.out, vals.out))
    # print( summary( pred_pm.dt[, 6:10]))
    write_fst( pred_pm.dt, fname_out)
    
    note <- paste( 'Unit conversions saved to', fname_out)
    message( note)
    return( note)
    
  } else {
    hyads.dt <- read.fst( fname, columns = c( 'x', 'y', 'uID', 'hyads'), as.data.table = T)
    hyads.dt <- hyads.dt[!is.na( x) & !is.na( y)]
    hyads.dt.c <- dcast( hyads.dt, x + y ~ uID, value.var = 'hyads')
    
    hyads.use.p <- rasterFromXYZ( hyads.dt.c, crs = p4s)
    hyads.use.p[is.na( hyads.use.p)] <- 0
    hyads.proj <- project_and_stack( dat.s[[1]], hyads.use.p, mask.use = mask.use)
    hyads.proj <- dropLayer( hyads.proj, 1)
    hyads.n <- paste0( 'X', names( hyads.dt.c)[!(names( hyads.dt.c) %in% c( 'x', 'y'))])
    names( hyads.proj) <- hyads.n
    
    # take fraction of total hyads
    hyads.frac.total <- hyads.proj / hyads.total.proj[[name.x]]
    hyads.frac.total.pm <- hyads.frac.total * pred_nomean.r
    names( hyads.frac.total.pm) <- hyads.n
    
    coords.out <- coordinates( hyads.frac.total.pm)
    vals.out <- values( hyads.frac.total.pm)
    
    # write out the data.table as fst
    pred_pm.dt <- data.table( cbind( coords.out, vals.out))
    # print( summary( pred_pm.dt[, 6:10]))
    write_fst( pred_pm.dt, fname_out)
    
    note <- paste( 'Unit conversions saved to', fname_out)
    message( note)
    return( note)
  }
}
