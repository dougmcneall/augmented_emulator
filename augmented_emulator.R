# augmented_emulator.R
# This code augments a Gaussian process emulator using temperature
# and precipitation (usually considered model outputs) as inputs.
# It allows biases in those fields to be accounted for when choosing
# plausible values of input parameters for the land surface component
# of the reduced-resolution climate model, FAMOUS.
# This is the main analysis code for:
# McNeall, D.J., Williams J., Betts, R.A. , Booth, B.B.B., Challenor, P.G. , Good, P.
# & Wiltshire A. (2019) Correcting a bias in a climate model with an augmented emulator,
# submitted as a discussion paper to Geoscientific Model Development
# Contact Doug McNeall dougmcneall@gmail.com @dougmcneall

# ------------------------------------------------------------
# Load packages and data
# ------------------------------------------------------------
library(DiceKriging)
library(RColorBrewer)
library(MASS)
library(fields)
library(parallel)
library(viridisLite)

load('famous_forest_fraction.RData')
load('famous_agg.Rdata')

# Load specific versions of a github repository
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/emtools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/imptools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/vistools.R')

# pallettes
rb <- brewer.pal(11, "RdBu")
ryg <- brewer.pal(11, "RdYlGn")
pbg <- brewer.pal(9, "PuBuGn")
bg <- brewer.pal(9, "BuGn")
yg <- brewer.pal(9, "YlGn")
byr <- rev(brewer.pal(11,'RdYlBu'))
br <- rev(rb)
blues <-  brewer.pal(9,'Blues')
rblues <-  rev(blues)

greens <-  brewer.pal(9,'Greens')
ygb <- brewer.pal(9, "YlGnBu")
brbg <- brewer.pal(11, "BrBG")
yob <- brewer.pal(9, "YlOrBr")
yor <- brewer.pal(9, "YlOrRd")
acc <- brewer.pal(8,'Paired')
cbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

col.amaz <- cbPal[1]
col.seasia <- cbPal[2]
col.congo <- cbPal[3]
col.global <- cbPal[4]
col.namerica <- cbPal[5]

col.tropics = c(rep(col.amaz, 100), rep(col.seasia,100), rep(col.congo,100))

pch.global <- 3
pch.amaz <- 1
pch.congo <- 2
pch.seasia <- 5
pch.namerica <- 4

# --------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------

dfunc.up <- function(x,y,...){
  # function for plotting 2d kernel density estimates in pairs() plot.
  require(MASS)
  require(RColorBrewer)
  
  br <- brewer.pal(9, 'Blues')
  kde <- kde2d(x,y)
  image(kde, col = br, add = TRUE)
}

dfunc.up.truth = function(x,y, ...){
  # function for plotting 2d kernel density estimates in pairs() plot,
  # adding a data point overlay.
  require(MASS)
  require(RColorBrewer)
  
  xtrue <- tail(x,1)
  ytrue <- tail(y,1)
  
  xdash <- head(x, -1)
  ydash <- head(y, -1)
  
  br <- brewer.pal(9, 'Blues')
  
  kde <- kde2d(xdash,ydash)
  image(kde, col = br, add = TRUE)
  points(xtrue, ytrue, pch =21, col = 'black', bg = 'red', cex = 1.5)
}

shadowtext <- function(x, y=NULL, labels, 
                       col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  # put text on a plot - white with a black background
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
          labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}

col3rd = function(n, pal, z){
  # produce a set of colours that match the values of z
  # Use for colour in 3rd dimension on a scatterplot.
  cRP  = colorRampPalette(pal)
  cols = cRP(n)
  
  #out = cols[(z - min(z))/diff(range(z))*n + 1 ]
  out = cols[cut(z, breaks = n) ]
  out
}

reset <- function() {
  # for putting labels outside the margins
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}


X = famous_agg[, 2:8]
X.norm = normalize(X)
X.stan.norm <- normalize(matrix(X.standard, nrow = 1), wrt = X)
colnames(X.stan.norm) = colnames(X)

# Build input matrix of (repeated) input parameters, modelled temperatures,
# precipitation and forest fractions.
X_AMAZON = famous_agg[ ,c(2,3,4,5,6,7,8,15,19)] # AMAZON
X_SEASIA = famous_agg[, c(2,3,4,5,6,7,8, 17,21)] # SEASIA
X_CONGO = famous_agg[, c(2,3,4,5,6,7,8,16,20)] # CONGO

colnames(X_AMAZON) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')
colnames(X_SEASIA) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')
colnames(X_CONGO) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')

X_tropics = rbind(X_AMAZON, X_SEASIA, X_CONGO)
X_tropics_norm = normalize(X_tropics)

# Output vector of broadleaf forest fraction
Y_tropics = c(famous_agg$AMAZ_MOD_FRAC, famous_agg$SEASIA_MOD_FRAC, famous_agg$CONGO_MOD_FRAC)

# Standard emulator fit to all tropical forests
tropics_fit = km(~., design = X_tropics_norm, response=Y_tropics)

# ----------------------------------------------------------------------------------
# Emulator diagnostics
#
# ----------------------------------------------------------------------------------
run_diagnostics = FALSE # The diagnostics section is slow, so only set to TRUE if you have the time
if(run_diagnostics){
  
  # Plot the emulator against the true leave-one-out prediction
  # This usually gives a very similar response to trend.reestim = TRUE
  # in leaveOneOut.km
  true.loo = function(X,y){
    out.mean = rep(NA, length(y))
    out.sd = rep(NA, length(y))
    
    for(i in 1:nrow(X)){
      X.trunc = X[-i, ]
      y.trunc = y[-i]
      
      X.target = matrix(X[i, ], nrow = 1)
      colnames(X.target) <- colnames(X)
      X.target.df = data.frame(X.target)
      
      fit = km(~., design = X.trunc, response = y.trunc)
      pred = predict(fit,newdata = X.target, type = 'UK')
      out.mean[i] = pred$mean
      out.sd[i] = pred$sd
    }
    return(list(mean = out.mean, sd = out.sd))
  }  
  
  true.loo.all = true.loo(X = X_tropics_norm, y = Y_tropics)
  
  true.loo.amazon = true.loo (X = X.norm, y = famous_agg$AMAZ_MOD_FRAC)
  true.loo.seasia = true.loo (X = X.norm, y = famous_agg$SEASIA_MOD_FRAC)
  true.loo.congo = true.loo (X = X.norm, y = famous_agg$CONGO_MOD_FRAC)
  
  mean(abs(true.loo.amazon$mean - famous_agg$AMAZ_MOD_FRAC))
  mean(abs(true.loo.seasia$mean - famous_agg$SEASIA_MOD_FRAC))
  mean(abs(true.loo.congo$mean - famous_agg$CONGO_MOD_FRAC))
  
  # Mean error when not using T/P emulator
  regular.mae = mean(abs(c(true.loo.amazon$mean, true.loo.seasia$mean, true.loo.congo$mean) - Y_tropics))
  print(paste('MAE for regular emulator =', regular.mae))
  
  aug.mae = mean(abs(true.loo.all$mean - Y_tropics))
  print(paste('MAE for augmented emulator =', aug.mae))
  
  # Mean absolute error is about 0.03 or 3%
  print(paste('With T/P mean absolute cross validation error = ', mean(abs(true.loo.all$mean - Y_tropics))))
  mean(abs(true.loo.all$mean[1:100] - famous_agg$AMAZ_MOD_FRAC))
  mean(abs(true.loo.all$mean[101:200] - famous_agg$SEASIA_MOD_FRAC))
  mean(abs(true.loo.all$mean[201:300] - famous_agg$CONGO_MOD_FRAC))
  
  pdf(width = 6, height = 6, file = 'graphics/true_loo_all.pdf' )
  xlim = c(-0.05, 1.05)
  ylim = c(-0.05, 1.05)
  par(las =1)
  plot(Y_tropics, true.loo.all$mean,
       pch = c(rep(21, 100), rep(22, 100), rep(24, 100)),
       bg = c(rep(col.amaz, 100), rep(col.seasia, 100), rep(col.congo, 100)),
       xlab = 'simulated forest fraction', ylab = 'emulated forest fraction',
       col = col.tropics,
       xlim = xlim, 
       ylim = ylim,
       bty = 'n',
       axes = FALSE,
       xaxs = 'i', yaxs = 'i')
  
  segments(x0 = Y_tropics, y0 = true.loo.all$mean - (2*true.loo.all$sd),
           x1 = Y_tropics, y1 = true.loo.all$mean +(2*true.loo.all$sd),
           col = col.tropics)
  axis(1, pos = 0, col = 'grey')
  axis(2, pos = 0, col = 'grey')
  abline(0,1, col = 'grey')
  legend('topleft', legend = c('Amazon', 'Asia', 'Africa'),
         inset = 0.1,
         pch = c(21, 22, 24), col = c(col.amaz, col.seasia, col.congo),
         pt.bg = c(col.amaz, col.seasia, col.congo),
         bty = 'n')
  
  legend('bottomright', legend = 'vertical lines depict \u00B1 2 \n standard deviations',
         inset = 0.1,
         pch = '',
         col = 'black',
         bty = 'n')
  dev.off()

  
  
  # Rank histograms for checking the uncertainty?
  # The principle behind the rank histogram is quite simple. 
  # Ideally, one property that is desired from an EF is reliable probabilities;
  # if ensemble relative frequency suggests P percent probability of occurrence,
  # the event truly ought to have P probability of occurring. 
  # For this probability to be reliable, the set of ensemble member forecast values
  # at a given point and the true state (the verification) ought to be able to be 
  # considered random samples from the same probability distribution.
  # This reliability then implies in turn that if an n-member ensemble and the
  # verification are pooled into a vector and sorted from lowest to highest,
  # then the verification is equally likely to occur in each of the n + 1 possible ranks. 
  # If the rank of the verification is tallied and the process repeated over many
  # independent sample points, a uniform histogram over the possible ranks should result.
  # From https://journals.ametsoc.org/doi/full/10.1175/1520-0493%282001%29129%3C0550%3AIORHFV%3E2.0.CO%3B2
  
  loo.rankhist = function(obs, pred.mean, pred.sd, n = 1000){
    # a version of the rank histogram     
    obs.ranks = rep(NA, length(obs))
    
    for(i in 1:length(obs)){
      ranks = rank(c(obs[i], rnorm(n = n, mean = pred.mean[i], sd = pred.sd[i])))
      obs.ranks[i] = ranks[1]
    }
    
    out = obs.ranks/(n+1)
    out
  }
  
  true.loo.rankhist = loo.rankhist(obs = Y_tropics, pred.mean = true.loo.all$mean, pred.sd = true.loo.all$sd, n = 500)
  
  pdf(file = 'graphics/rankhist.pdf', width = 6, height =4)
  par(las = 1, fg = 'white')
  hist(true.loo.rankhist, col = 'grey', 
       axes = FALSE, main = '',
       xlab = 'Rank'
  )
  axis(1, col = 'black')
  axis(2, col = 'black')
  dev.off()
  
  # Prediction errors normalised by their standard deviations
  # approximately follow a normal distribution - pehaps the tails are a little off.
  true.loo.err.norm = (true.loo.all$mean - Y_tropics) / true.loo.all$sd
  pdf(file = 'graphics/normalisedQQ.pdf')
  par(las = 1)
  qqnorm(true.loo.err.norm)
  abline(0,1)
  dev.off()
  
}else{print('skipping diagnostics')}

# ------------------------------------------------------------------
# Bias correction section
# ------------------------------------------------------------------

# normalize amazon obs
tp.amaz.norm <- normalize(
  matrix(c(temps_obs$AMAZ_OBS_TEMP+273.15, precips_obs$AMAZ_OBS_PRECIP),nrow=1),
  wrt=X_tropics[, 8:9]
)
#colnames(tp.amaz.norm) <- c('OBS_TEMP', 'OBS_PRECIP')
model_climate_names = c('MOD_TEMP', 'MOD_PRECIP')
colnames(tp.amaz.norm) <- model_climate_names

tp.seasia.norm <- normalize(
  matrix(c(temps_obs$SEASIA_OBS_TEMP+273.15, precips_obs$SEASIA_OBS_PRECIP),nrow=1),
  wrt=X_tropics[, 8:9]
)
#colnames(tp.seasia.norm) <- c('OBS_TEMP', 'OBS_PRECIP')
colnames(tp.seasia.norm) <- model_climate_names

tp.congo.norm <- normalize(
  matrix(c(temps_obs$CONGO_OBS_TEMP+273.15, precips_obs$CONGO_OBS_PRECIP),nrow=1),
  wrt=X_tropics[, 8:9]
)
#colnames(tp.congo.norm) <- c('OBS_TEMP', 'OBS_PRECIP')
colnames(tp.congo.norm) <- model_climate_names

# Default parameters attached and observed temperature and precip
amaz.x  <- cbind(X.stan.norm, tp.amaz.norm)
congo.x <- cbind(X.stan.norm, tp.congo.norm)
seasia.x <- cbind(X.stan.norm, tp.seasia.norm)

# Emulator predicted bias corrected at default parameters
pred.amaz.bc <- predict(tropics_fit, newdata=amaz.x, type='UK')
pred.congo.bc <- predict(tropics_fit, newdata=congo.x, type='UK')
pred.seasia.bc <- predict(tropics_fit, newdata=seasia.x, type='UK')

# Emulator fit with no temp/precip
fit.amazon <- km(~., design=X.norm, response=famous_agg$AMAZ_MOD_FRAC)
fit.congo <- km(~., design=X.norm, response=famous_agg$CONGO_MOD_FRAC)
fit.seasia <- km(~., design=X.norm, response=famous_agg$SEASIA_MOD_FRAC)

standard.amazon <- predict(fit.amazon, newdata=X.stan.norm, type='UK')
standard.congo <- predict(fit.congo, newdata=X.stan.norm, type='UK')
standard.seasia <- predict(fit.seasia, newdata=X.stan.norm, type='UK')

obs_amazon = obs[,'AMAZON']
obs_congo = obs[,'CONGO']
obs_seasia = obs[,'SEASIA']

# It looks like bias correction via T and P improves the Amazon,
# slightly worsens the Congo, and very slightly improves SE Asia. 
obs_amazon - standard.amazon$mean
obs_amazon - pred.amaz.bc$mean

obs_congo - standard.congo$mean
obs_congo - pred.congo.bc$mean

obs_seasia - standard.seasia$mean
obs_seasia- pred.seasia.bc$mean

# --------------------------------------------------------------------
# Sensitivity analysis, including temperature and precipitation
# How does the forest fraction sensitivity to parameters change
# at the default settings for all climates?
# --------------------------------------------------------------------
library(sensitivity)

n = 21
xlist = list(amaz.x, seasia.x, congo.x)
ylist = list(famous_agg$AMAZ_MOD_FRAC, famous_agg$SEASIA_MOD_FRAC, famous_agg$CONGO_MOD_FRAC)
## build a matrix of OAT predictions
oat.mean.mat = matrix(nrow = n*length(amaz.x), ncol = length(xlist))
oat.sd.mat = matrix(nrow = n*length(amaz.x), ncol = length(xlist))
fit.sens = km(~., design = X_tropics_norm, response = Y_tropics)

for(i in 1:length(xlist)){
  
  X.oat = oaat.design(X_tropics_norm, n = n, hold = xlist[[i]])
  colnames(X.oat) = colnames(xlist[[i]])
  pred.sens = predict(fit.sens, newdata = X.oat, type = 'UK')
  oat.mean.mat[, i ] = pred.sens$mean
  oat.sd.mat[, i ] = pred.sens$sd
}

col.list = list(col.amaz, col.seasia, col.congo)

pdf(width = 7, height = 6, file = 'graphics/sensitivity_TP_all.pdf') 
par(mfrow = c(2,5), las = 1, mar = c(5,0.5,3,0.5), oma = c(0,5,0,0), fg = 'grey')

for(i in 1: ncol(X_tropics_norm)){
  
  ix <- seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oat[ix, i], pred.sens$mean[ix], ylim = c(0,1), xlab = colnames(X.oat)[i], type = 'n', axes = FALSE)
  axis(1)
  if (i==1 | i==6 ) {axis(2)
    mtext(side = 2, line = 3.5, text = 'Forest fraction', las = 0, col = 'black')
  }
  for(j in 1:length(xlist)){
    
    col.chosen = col.list[[j]]
    col.transp = adjustcolor(col.chosen, alpha = 0.5)
    
    pred.mean = oat.mean.mat[, j]
    pred.sd   = oat.sd.mat[, j]
    
    polygon(x = c(X.oat[ix, i], rev(X.oat[ix, i])),
            y = c( (pred.mean[ix] - pred.sd[ix]), rev(pred.mean[ix] + pred.sd[ix])),
            col = col.transp, border = col.transp)
    
    lines(X.oat[ix, i], pred.mean[ix], ylim = c(0,1), xlab = colnames(X.oat)[i], col = col.chosen )
  }
}
reset()
legend('top', legend = c('Amazon', 'SE Asia', 'C Africa'), 
       col = c(col.list, recursive = TRUE),
       lty = 'solid', lwd = 1, pch = NA, bty = 'n',
       text.col = 'black',
       fill = adjustcolor(c(col.list, recursive = TRUE), alpha = 0.5),
       cex = 1.2, border = NA, horiz = TRUE)

dev.off()

# ------------------------------------------------------------------
# FAST99 sensitivity analysis of Saltelli et al (1999)
# generate the design to run the emulator at, using fast99
# ------------------------------------------------------------------

X.fast <- fast99(model = NULL, factors = colnames(X_tropics_norm), n = 1000,
                 q = "qunif", q.arg = list(min = 0, max = 1))

pred.fast = predict(tropics_fit, newdata = X.fast$X, type = 'UK')



# Calculate the sensitivity indices
fast.tell <- tell(X.fast, pred.fast$mean)

bp.convert <- function(fastmodel){
  # get the FAST summary into an easier format for barplot
  fast.summ <- print(fastmodel)
  fast.diff <- fast.summ[ ,2] - fast.summ[ ,1]
  fast.bp <- t(cbind(fast.summ[ ,1], fast.diff))
  fast.bp
}

pdf(width = 7, height = 5, file = 'graphics/fast_barplot.pdf')
par(las = 2, mar = c(9,5,3,2))
barplot(bp.convert(fast.tell), col = c('skyblue', 'grey'), ylab = 'relative sensitivity')
legend('topleft',legend = c('Main effect', 'Interactions'), fill = c('skyblue', 'grey') )
dev.off()


# Reviewer 1 asks how the Fast algorithm might be impacted because some of the
# FAST99 design points are in a zone with few design points (eg, the cool, wet
# corner of T/P space.

# First, what is the implausibility of the FAST99 X?

fast.impl.amaz = impl(em = pred.fast$mean, em.sd = pred.fast$sd,
                  disc = 0, obs = obs_amazon, disc.sd = 0.01, obs.sd = 0)


fast.impl.seasia= impl(em = pred.fast$mean, em.sd = pred.fast$sd,
                  disc = 0, obs = obs_seasia, disc.sd = 0.01, obs.sd = 0)


fast.impl.congo = impl(em = pred.fast$mean, em.sd = pred.fast$sd,
                  disc = 0, obs = obs_congo, disc.sd = 0.01, obs.sd = 0)

dev.new()
par(mfrow = c(3,1))
xlim = c(0,12)
hist(fast.impl.amaz, xlim = xlim)
hist(fast.impl.seasia, xlim = xlim)
hist(fast.impl.congo, xlim = xlim)

dev.new()
plot(X.fast$X[,'MOD_TEMP'], X.fast$X[,'MOD_PRECIP'])

# Calculating the 'Shapley measures', which allows for non-independence
# of the inputs

# This doesn't give great results
#shapMC = shapleySubsetMc(X=X_tropics,Y=Y_tropics,
#  Ntot=NULL, Ni=3, cat=NULL, weight=NULL, discrete=NULL)

#fast.firstorder = print(fast.tell)[,1]
#fast.total = print(fast.tell)[,2]

#dev.new()
#plot(fast.firstorder, shapMC$shapley)


#X.fast <- fast99(model = NULL, factors = colnames(X_tropics_norm), n = 1000,
#                 q = "qunif", q.arg = list(min = 0, max = 1))

#pred.fast = predict(tropics_fit, newdata = X.fast$X, type = 'UK')




#shapleyPermRand(model = NULL, Xall, Xset, d, Nv, m, No = 1, Ni = 3, colnames = NULL, ...)

     ## S3 method for class 'shapleyPermRand'

#X.shapPR = shapleyPermRand(model = NULL, Xall, Xset, d, Nv, m, No = 1, Ni = 3, colnames = NULL, ...)

#pred.shapPR = predict(tropics_fit, newdata = X.fast$X, type = 'UK')

#tell(x, y = NULL, return.var = NULL, ...)



# ------------------------------------------------------
# Find the set of plausible inputs, when 
# temperature and precip are included in the inputs
# ------------------------------------------------------

inputs.set <- function(X, y, thres, obs, obs.sd = 0, disc = 0, disc.sd = 0, n = 100000, abt = FALSE){ 
  # find a set of inputs that are consistent with a particular
  # set of implausibility (either below or above)
  
  X.mins <- apply(X,2,min)
  X.maxes <- apply(X,2,max)
  X.unif <- samp.unif(n, mins = X.mins, maxes = X.maxes)
  colnames(X.unif) <- colnames(X)
  
  fit <- km(~., design = X, response = y, control = list(trace = FALSE))
  pred <- predict(fit, newdata = X.unif, type = 'UK')
  pred.impl <- impl(em = pred$mean, em.sd = pred$sd,
                    disc = disc, obs = obs, disc.sd = disc.sd, obs.sd = obs.sd)
  
  if(abt){
    # choose those above the threshold 
    ix.bt <- pred.impl > thres
  }
  
  else{
    ix.bt <- pred.impl < thres
  }
  
  X.out <- X.unif[ix.bt, ]
  
  return(list(X.out = X.out, fit = fit, X.unif = X.unif, pred = pred,pred.impl = pred.impl))   
}


plausible.amazon.bc <- inputs.set(X = X_tropics_norm, y = Y_tropics,thres = 3,
                                  obs = obs_amazon,
                                  obs.sd = 0,
                                  disc = 0,
                                  disc.sd = 0.01,
                                  n = 100000,
                                  abt = FALSE)

# Emulated surface of forest fraction given T & P
#dev.new(width = 5, height =5)
pdf(file = 'graphics/emulated_fraction_vs_temp_precip.pdf',width = 7, height = 7)
quilt.plot(plausible.amazon.bc$X.unif[,8], plausible.amazon.bc$X.unif[, 9], plausible.amazon.bc$pred$mean, 
           col = byr, xlab = 'Normalised regional temperature', ylab = 'Normalised regional precipitation',
           legend.args = list(text = "forest\nfraction",col="black", cex=1.2, side=3, line=1))
cex = 1.5
points(X_tropics_norm[1:100,8], X_tropics_norm[1:100,9], col = 'black', bg = col.amaz, pch = 21, cex = cex)
points(X_tropics_norm[101:200,8], X_tropics_norm[101:200,9], col = 'black', bg = col.seasia, pch = 21, cex = cex)
points(X_tropics_norm[201:300,8], X_tropics_norm[201:300,9], col = 'black', bg = col.congo, pch = 21, cex = cex)

points(tp.amaz.norm, col = 'black', pch = 21, cex = 2.5, bg = col.amaz, lwd = 2.5)
points(tp.congo.norm, col = 'black', pch = 21, cex = 2.5, bg = col.congo, lwd = 2.5)
points(tp.seasia.norm, col = 'black', pch = 21, cex = 2.5, bg = col.seasia, lwd = 2.5)

text(tp.amaz.norm, 'Amazon', pos = 4, font = 2)
text(tp.congo.norm, 'Central Africa', pos = 4, font = 2)
text(tp.seasia.norm, 'SE Asia', pos = 4, font = 2)
dev.off()

# No emulated surface in this version
pdf(file = 'graphics/fraction_vs_temp_precip_pcolcor.pdf',width = 8, height = 7)
par(las = 1, fg = 'black', mar = c(5,6,3,7))
Y_obs = c(Y_tropics, obs_amazon, obs_seasia, obs_congo)
Y_obs_ix = 1:300
zcolor = col3rd(n=9, pal= viridis(9), z = Y_obs) 
pr = 1e5 # precipitation scaling factor
plot(X_tropics[, 8]-273.15, X_tropics[, 9]*pr, col = 'black', pch = 21, cex = 2, type = 'n',
     xlab = expression(paste('Regional temperature (', degree,'C)')),
     ylab = expression(paste('Regional Precipitation x',10^5,' kgm'^-2,'s'^-1))
)

cex = 1.4
lwd = 1.5

points(X_tropics[1:100,8]-273.15, X_tropics[1:100,9]*pr, col = 'black', bg = zcolor[1:100], pch = 21, cex = cex, lwd = lwd)
points(X_tropics[101:200,8]-273.15, X_tropics[101:200,9]*pr, col = 'black', bg = zcolor[101:200], pch = 22, cex = cex, lwd = lwd)
points(X_tropics[201:300,8]-273.15, X_tropics[201:300,9]*pr, col = 'black', bg = zcolor[201:300], pch = 24, cex = cex, lwd = lwd)

points(temps_obs$AMAZ_OBS_TEMP, precips_obs$AMAZ_OBS_PRECIP*pr, col = 'black', pch = 21, cex = 2.5, bg = zcolor[301], lwd = 2)
points(temps_obs$SEASIA_OBS_TEMP, precips_obs$SEASIA_OBS_PRECIP*pr, col = 'black', pch = 22, cex = 2.5, bg = zcolor[302], lwd = 2)
points(temps_obs$CONGO_OBS_TEMP, precips_obs$CONGO_OBS_PRECIP*pr, col = 'black', pch = 24, cex = 2.5, bg = zcolor[303], lwd = 2)

shadowtext(temps_obs$AMAZ_OBS_TEMP,precips_obs$AMAZ_OBS_PRECIP*pr, 'Amazon', pos = 4, font = 2,r =0.2)
shadowtext(temps_obs$SEASIA_OBS_TEMP, precips_obs$SEASIA_OBS_PRECIP*pr, 'SE Asia', pos = 4, font = 2, r = 0.2)
shadowtext(temps_obs$CONGO_OBS_TEMP, precips_obs$CONGO_OBS_PRECIP*pr, 'Central Africa', pos = 4, font = 2, r = 0.2)
legend('topleft', pch = c(21,22,24,24,24), pt.bg = c(NA,NA,NA,NA,viridis(9)[9]), legend = c('Amazon', 'SE Asia', 'Central Africa', 'large points are observations', 'fill colour is forest fraction'),
       pt.cex = c(1,1,1,1.4,1.4), bty = 'n')

image.plot(z = Y_obs, legend.only = TRUE, col = viridis(9), horizontal = FALSE,  legend.args = list(text = "forest\nfraction",col="black", cex=1.2, side=3, line=1))
dev.off()

# ------------------------------------------------------------------------
# Two-at-a-time sensitivity analysis, with inputs held at their default
# settings.
# ------------------------------------------------------------------------

# construct output matrix
tseq = seq(from = 0, to = 1, by = 0.05)
pseq = seq(from = 0, to = 1, by = 0.05)
tdashp = expand.grid(tseq, pseq)
norm.mat = matrix(X.stan.norm, nrow = nrow(tdashp), ncol = ncol(X.stan.norm), byrow = TRUE)
X.taat.tp = cbind(norm.mat, tdashp)
colnames(X.taat.tp) = colnames(X_tropics_norm)

# sample from the emulator
y.taat.tp = predict(tropics_fit, newdata = X.taat.tp, type = 'UK')

# Can use quilt plot as long as the number of points in either direction matches the
# data.
ncol = 11
allz = c(Y_tropics,obs_amazon,obs_seasia, obs_congo, y.taat.tp$mean)
zcolor = col3rd(n=ncol, pal=viridis(ncol), z = allz) 

pdf(width = 7, height = 7, file = 'graphics/taat_temp_precip_quilt.pdf')
par(las = 1)
quilt.plot(X.taat.tp[,8], X.taat.tp[,9], y.taat.tp$mean, 
           col = viridis(ncol), nx = 21, ny = 21,
           xlab = 'Normalised Regional Mean Temperature', ylab = 'Normalised Regional Mean Precipitation')

cex = 1.4
lwd = 1.5
points(X_tropics_norm[1:100,8], X_tropics_norm[1:100,9], 
       col = 'black', bg = zcolor[1:100], pch = 21, cex = cex, lwd = lwd)
points(X_tropics_norm[101:200,8], X_tropics_norm[101:200,9], 
       col = 'black', bg = zcolor[101:200], pch = 22, cex = cex, lwd = lwd)
points(X_tropics_norm[201:300,8], X_tropics_norm[201:300,9], 
       col = 'black', bg = zcolor[201:300], pch = 24, cex = cex, lwd = lwd)

points(tp.amaz.norm, col = 'black', pch = 21, cex = 2.5, bg = zcolor[301], lwd = 2)
points(tp.seasia.norm, col = 'black', pch = 22, cex = 2.5, bg = zcolor[302], lwd = 2)
points(tp.congo.norm, col = 'black', pch = 24, cex = 2.5, bg = zcolor[303], lwd = 2)

shadowtext(tp.amaz.norm[1],tp.amaz.norm[2], 'Amazon', pos = 4, font = 2,r =0.2)
shadowtext(tp.congo.norm[1], tp.congo.norm[2], 'Central Africa', pos = 4, font = 2, r = 0.2)
shadowtext(tp.seasia.norm[1], tp.seasia.norm[2], 'SE Asia', pos = 4, font = 2, r = 0.2)

dev.off()

# ------------------------------------------------------------------------
# We worked out the response surface for y = f(X, T, P).
# Now bias correct T and P in all of the ensemble members.
# ------------------------------------------------------------------------

# Fit the whole data set
fit.tropics = km(~.,design = X_tropics_norm, response = Y_tropics)

# Bias correct the model runs using the observed (normalised)
# temperature and precipitation.
bc.x  = rbind(cbind(X.norm, matrix(tp.amaz.norm, nrow = 100, ncol = 2, byrow = TRUE)),
              cbind(X.norm, matrix(tp.seasia.norm, nrow = 100, ncol = 2, byrow = TRUE)),
              cbind(X.norm, matrix(tp.congo.norm, nrow = 100, ncol = 2, byrow = TRUE))
)
colnames(bc.x) <- colnames(X_tropics_norm)

model_runs.bc = predict(fit.tropics, newdata = bc.x , type = 'UK')

plot(Y_tropics,model_runs.bc$mean, xlim = c(0,1), ylim = c(0,1))
abline(0,1)

dev.new()
par(mfrow = c(2,3))
hist(Y_tropics[1:100] - obs_amazon, xlim = c(-1,1))
hist(Y_tropics[101:200] - obs_seasia, xlim = c(-1,1))
hist(Y_tropics[201:300] - obs_congo, xlim = c(-1,1))

hist(model_runs.bc$mean[1:100] - obs_amazon, xlim = c(-1,1))
hist(model_runs.bc$mean[101:200] - obs_seasia, xlim = c(-1,1))
hist(model_runs.bc$mean[201:300] - obs_congo, xlim = c(-1,1))

# How much does the forest fraction improve?
mean(abs(Y_tropics[1:100] - obs_amazon))
mean(abs(model_runs.bc$mean[1:100] - obs_amazon))

mean(abs(Y_tropics[101:200] - obs_seasia))
mean(abs(model_runs.bc$mean[101:200] - obs_seasia))

mean(abs(Y_tropics[201:300] - obs_congo))
mean(abs(model_runs.bc$mean[201:300] - obs_congo))

# Plot the observation forest fraction, and the model at the default 
# parameters and bias corrected.
pdf(width = 7, height = 5, file='graphics/bias_corrected_fractions.pdf')
par(las=1)
plot(c(1,2,3), c(obs_amazon, obs_congo, obs_seasia), xlim=c(0.5,3.5), ylim=c(0,1), pch=19,
     col = c(col.amaz, col.congo, col.seasia), cex=1.5, axes=FALSE, xlab='', ylab='Forest fraction',
     pty = 'n')

points(rep(1.1, 100), Y_tropics[1:100], col = 'grey', pch = '-')
points(rep(1.2, 100), model_runs.bc$mean[1:100], col = 'pink', pch = '-')

points(rep(2.1, 100), Y_tropics[101:200], col = 'grey', pch = '-')
points(rep(2.2, 100), model_runs.bc$mean[101:200], col = 'pink', pch = '-')

points(rep(3.1, 100), Y_tropics[201:300], col = 'grey', pch = '-')
points(rep(3.2, 100), model_runs.bc$mean[201:300], col = 'pink', pch = '-')

points(c(1,2,3), c(obs_amazon, obs_congo, obs_seasia), pch=19,
       col = c(col.amaz, col.congo, col.seasia), cex=1.5)

points(c(1.1,2.1,3.1), c(standard.amazon$mean, standard.congo$mean,standard.seasia$mean), pch=19)

segments(1.1, standard.amazon$mean - standard.amazon$sd,
         1.1, standard.amazon$mean + standard.amazon$sd, col='black')

segments(2.1, standard.congo$mean - standard.congo$sd,
         2.1, standard.congo$mean + standard.congo$sd, col='black')

segments(3.1, standard.seasia$mean - standard.seasia$sd,
         3.1, standard.seasia$mean + standard.seasia$sd, col='black')

points(c(1.2,2.2,3.2),c(pred.amaz.bc$mean, pred.congo.bc$mean, pred.seasia.bc$mean),  col='red3', pch=19)

segments(1.2, pred.amaz.bc$mean - pred.amaz.bc$sd,
         1.2, pred.amaz.bc$mean + standard.amazon$sd, col='red3')

segments(2.2, pred.congo.bc$mean - standard.congo$sd,
         2.2, pred.congo.bc$mean + standard.congo$sd, col='red3')

segments(3.2, pred.seasia.bc$mean - pred.seasia.bc$sd,
         3.2, pred.seasia.bc$mean + pred.seasia.bc$sd, col='red3')

axis(1, labels = c('Amazon', 'Africa', 'SE Asia'), at = c(1.1,2.1,3.1))
axis(2)

text(1, obs_amazon, 'observation', pos=2, col='grey', cex=0.9)
text(1.1, standard.amazon$mean, 'default\nparameters', pos=2, cex=0.9, col='grey')
text(1.2, pred.amaz.bc$mean, 'bias\ncorrected', pos=4, col='red3', cex=0.9)
legend('topright', legend = c('model runs', 'bias corrected model runs'),
       col = c('grey', 'pink'), pch = '-', bty = 'n'
)
dev.off()


pdf(file = 'graphics/dotchart_fractions.pdf', width = 6, height = 7)
par(mar = c(5,6,2,8))
plot(c(obs_amazon, obs_congo, obs_seasia), c(1,2,3), ylim=c(0.5,3.5), xlim=c(0,1), pch=19,
     col = 'blue', cex=1.5, ylab='', xlab='Forest fraction',
     axes = FALSE, type = 'n', bty = 'l',
     panel.first = rect(c(0,0),c(0.5,2.5),c(1,1),c(1.5,3.5), col = 'grey95', border = 'grey95'),
     xaxs = 'i', yaxs = 'i'
)
abline(h = c(0.75, 1, 1.25, 1.75, 2, 2.25, 2.75, 3, 3.25), col = 'grey', lty = 'dashed')

# observed points
points(c(obs_amazon, obs_congo, obs_seasia),c(1.25,2.25,3.25), pch=19,
       col = c(col.amaz, col.congo, col.seasia), cex=1.5)

# standard model run points
points(c(standard.amazon$mean, standard.congo$mean,standard.seasia$mean),
       1:3, pch=17, col = c(col.amaz, col.congo, col.seasia), cex = 1.5)

segments(standard.amazon$mean - standard.amazon$sd,1,
         standard.amazon$mean + standard.amazon$sd,1, col=col.amaz, lwd = 1.5)

segments(standard.congo$mean - standard.congo$sd,2,
         standard.congo$mean + standard.congo$sd,2, col=col.congo, lwd =1.5)

segments(standard.seasia$mean - standard.seasia$sd,3,
         standard.seasia$mean + standard.seasia$sd,3, col=col.seasia, lwd = 1.5)

# bias corrected points
points(c(pred.amaz.bc$mean, pred.congo.bc$mean, pred.seasia.bc$mean),
       c(0.75,1.75,2.75), col=c(col.amaz, col.congo, col.seasia), pch=15, cex = 1.5)

segments( pred.amaz.bc$mean - pred.amaz.bc$sd, 0.75,
          pred.amaz.bc$mean + pred.amaz.bc$sd, 0.75, col=col.amaz, lwd = 1.5)

segments(pred.congo.bc$mean - pred.congo.bc$sd,1.75,
         pred.congo.bc$mean + pred.congo.bc$sd,1.75,  col=col.congo, lwd = 1.5)

segments(pred.seasia.bc$mean - pred.seasia.bc$sd,2.75,
         pred.seasia.bc$mean + pred.seasia.bc$sd,2.75, col=col.seasia, lwd = 1.5)

axis(2, labels = c('Amazon', 'Africa', 'SE Asia'), cex.axis = 1.3, at = 1:3, col = NA, las = 1)
axis(1)
par(xpd = TRUE)
text(c(1,1,1), c(3.25, 3, 2.75), labels = c('observed', 'default parameters', 'bias corrected'), pos = 4, cex = 1)
dev.off()

# Error for the standard runs when bias corrected and not
print(paste('amazon default error =', standard.amazon$mean - obs_amazon))
print(paste('SE Asia default error =', standard.seasia$mean - obs_seasia))
print(paste('C Africa default error =', standard.congo$mean - obs_congo))

print(paste('amazon bias corrected error =', pred.amaz.bc$mean - obs_amazon))
print(paste('SE Asia bias corrected error =', pred.seasia.bc$mean - obs_seasia))
print(paste('C Africa bias corrected error =', pred.congo.bc$mean - obs_congo))

nobc.mae = mean(abs(c((standard.amazon$mean - obs_amazon), (standard.seasia$mean - obs_seasia), (standard.congo$mean - obs_congo)))) 
bc.mae = mean(abs(c((pred.amaz.bc$mean - obs_amazon), (pred.seasia.bc$mean - obs_seasia), (pred.congo.bc$mean - obs_congo))))

bc.mae / nobc.mae

# --------------------------------------------------------
# Find points which are NROY for all three systems,
# When T and P are held at observed values and not 
# --------------------------------------------------------

# Find the set of plausible inputs, when temperature and precip are included in the inputs
# create a 'core' X
n = 100000
X.mins = apply(X.norm,2,min)
X.maxes = apply(X.norm,2,max)
X.unif = samp.unif(n, mins = X.mins, maxes = X.maxes)

# append the observed T and P
X.unif.amaz = cbind(X.unif, 
                    matrix(tp.amaz.norm, nrow = n,ncol = 2, byrow = TRUE))
colnames(X.unif.amaz) <- colnames(X_tropics_norm)
X.unif.seasia = cbind(X.unif, 
                      matrix(tp.seasia.norm, nrow = n,ncol = 2, byrow = TRUE))
colnames(X.unif.seasia) <- colnames(X_tropics_norm)
X.unif.congo = cbind(X.unif, 
                     matrix(tp.congo.norm, nrow = n,ncol = 2, byrow = TRUE))
colnames(X.unif.congo) <- colnames(X_tropics_norm)

disc = 0
obs.sd = 0
disc.sd = 0.01
thres = 3

pred.unif.amaz = predict(fit.amazon, newdata = X.unif, type = 'UK')
amaz.impl = impl(em = pred.unif.amaz$mean, em.sd = pred.unif.amaz$sd,
                 disc = disc, obs = obs_amazon,
                 disc.sd = disc.sd,
                 obs.sd = obs.sd)
nroy.ix.amaz = which(amaz.impl < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_amazon.pdf')
pairs(rbind(X.unif[nroy.ix.amaz, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

pred.unif.amaz.bc = predict(fit.tropics, newdata = X.unif.amaz, type = 'UK')
amaz.impl.bc = impl(em = pred.unif.amaz.bc$mean, em.sd = pred.unif.amaz.bc$sd,
                    disc = disc, obs = obs_amazon,
                    disc.sd = disc.sd,
                    obs.sd = obs.sd)
nroy.ix.amaz.bc = which(amaz.impl.bc < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_amazon_bc.pdf')
pairs(rbind(X.unif[nroy.ix.amaz.bc, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

# Now SE Asia
pred.unif.seasia = predict(fit.seasia, newdata = X.unif, type = 'UK')
seasia.impl = impl(em = pred.unif.seasia$mean, em.sd = pred.unif.seasia$sd,
                   disc = disc, obs = obs_seasia,
                   disc.sd = disc.sd,
                   obs.sd = obs.sd)
nroy.ix.seasia = which(seasia.impl < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_seasia.pdf')
pairs(rbind(X.unif[nroy.ix.seasia, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

pred.unif.seasia.bc = predict(fit.tropics, newdata = X.unif.seasia, type = 'UK')
seasia.impl.bc = impl(em = pred.unif.seasia.bc$mean, em.sd = pred.unif.seasia.bc$sd,
                      disc = disc, obs = obs_seasia,
                      disc.sd = disc.sd,
                      obs.sd = obs.sd)
nroy.ix.seasia.bc = which(seasia.impl.bc < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_seasia_bc.pdf')
pairs(rbind(X.unif[nroy.ix.seasia.bc, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

# Now Congo
pred.unif.congo = predict(fit.congo, newdata = X.unif, type = 'UK')
congo.impl = impl(em = pred.unif.congo$mean, em.sd = pred.unif.congo$sd,
                  disc = disc, obs = obs_congo,
                  disc.sd = disc.sd,
                  obs.sd = obs.sd)
nroy.ix.congo = which(congo.impl < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_congo.pdf')
pairs(rbind(X.unif[nroy.ix.congo, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

pred.unif.congo.bc = predict(fit.tropics, newdata = X.unif.congo, type = 'UK')
congo.impl.bc = impl(em = pred.unif.congo.bc$mean, em.sd = pred.unif.congo.bc$sd,
                     disc = disc, obs = obs_congo,
                     disc.sd = disc.sd,
                     obs.sd = obs.sd)
nroy.ix.congo.bc = which(congo.impl.bc < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_congo_bc.pdf')
pairs(rbind(X.unif[nroy.ix.congo.bc, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

# Implausibility at the default settings
amaz.impl.default = impl(em = standard.amazon$mean, em.sd = standard.amazon$sd,
                         disc = disc, obs = obs_amazon,
                         disc.sd = disc.sd,
                         obs.sd = obs.sd)

seasia.impl.default = impl(em = standard.seasia$mean, em.sd = standard.seasia$sd,
                           disc = disc, obs = obs_seasia,
                           disc.sd = disc.sd,
                           obs.sd = obs.sd)

congo.impl.default = impl(em = standard.congo$mean, em.sd = standard.congo$sd,
                          disc = disc, obs = obs_congo,
                          disc.sd = disc.sd,
                          obs.sd = obs.sd)

# Bias corrected Implausibility
amaz.impl.bc.default = impl(em = pred.amaz.bc$mean, em.sd = pred.amaz.bc$sd,
                            disc = disc, obs = obs_amazon,
                            disc.sd = disc.sd,
                            obs.sd = obs.sd)

seasia.impl.bc.default = impl(em = pred.seasia.bc$mean, em.sd = pred.seasia.bc$sd,
                              disc = disc, obs = obs_seasia,
                              disc.sd = disc.sd,
                              obs.sd = obs.sd)

congo.impl.bc.default = impl(em = pred.congo.bc$mean, em.sd = pred.congo.bc$sd,
                             disc = disc, obs = obs_congo,
                             disc.sd = disc.sd,
                             obs.sd = obs.sd)

# ----------------------------------------------------------------
# What part of parameter space matches everything?
# (We use a very low uncertainty)
# ----------------------------------------------------------------

nroy.bc.ix = intersect(intersect(nroy.ix.amaz.bc,nroy.ix.seasia.bc ), nroy.ix.congo.bc)
nroy.nobc.ix = intersect(intersect(nroy.ix.amaz,nroy.ix.seasia ), nroy.ix.congo)

(length(nroy.nobc.ix)/n) *100 # 2% of space is NROY with no bias correction.
(length(nroy.bc.ix)/n) *100 # 28% of space is NROY with bias correction.

pdf(width = 7, height = 7, file = 'graphics/best_inputs_all_bc.pdf')
pairs(rbind(X.unif[nroy.bc.ix, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

pdf(width = 7, height = 7, file = 'graphics/best_inputs_all_nobc.pdf')
pairs(rbind(X.unif[nroy.nobc.ix, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

# NROY inputs in the bias corrected case 
X.nroy.bc = rbind(X.unif[nroy.bc.ix, ], X.stan.norm)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_all_bc_default.pdf')
pairs(X.nroy.bc, panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

# Proportion of NROY space that is "shared"
prop.shared = function(a,b){
  out = length(intersect(a,b)) / length(union(a,b))
  out
}

a = nroy.ix.amaz.bc
b = nroy.ix.seasia.bc
d = nroy.ix.congo.bc

# bias corrected "proportion of shared space"
length(intersect(d,intersect(a,b))) / length(union(d,union(a,b)))

# as a proportion of total space:
length(intersect(d,intersect(a,b))) / n 

# non bias corrected "proportion of shared space"
a = nroy.ix.amaz
b = nroy.ix.seasia
d = nroy.ix.congo

length(intersect(d,intersect(a,b))) / length(union(d,union(a,b)))
# as a proportion of total space:
length(intersect(d,intersect(a,b))) / n 

# shared space with each pair of forests
# Bias-corrected
prop.shared(nroy.ix.amaz.bc, nroy.ix.seasia.bc)
prop.shared(nroy.ix.amaz.bc, nroy.ix.congo.bc)
prop.shared(nroy.ix.seasia.bc, nroy.ix.congo.bc)

# Non bias-corrected
prop.shared(nroy.ix.amaz, nroy.ix.seasia)
prop.shared(nroy.ix.amaz, nroy.ix.congo)
prop.shared(nroy.ix.seasia, nroy.ix.congo)

# Multicriteria optimisation?
# what does the parameter space with the lowest absolute error look like?
# the 'best' region?
hist(pred.unif.congo.bc$mean - obs_congo)
hist(pred.unif.seasia.bc$mean - obs_seasia)
hist(pred.unif.amaz.bc$mean - obs_amazon)

# Where is the error smallest? 

congo.ae = abs(pred.unif.congo.bc$mean - obs_congo)
seasia.ae = abs(pred.unif.seasia.bc$mean - obs_seasia)
amaz.ae = abs(pred.unif.amaz.bc$mean - obs_seasia)

total.ae = congo.ae + seasia.ae + amaz.ae

best.ix = which(total.ae < 0.25)

pdf(width = 7, height = 7, file = 'graphics/smallest_ae_inputs_all_bc_default.pdf')
pairs(rbind(X.unif[best.ix, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

# what does the output at those points look like?
# Even with a bias correction, the emulator indicates that 
# we underestimate congo and amazon, and overestimate SE Asia

pdf(width = 7, height = 7, file = 'graphics/smallest_ae_hists.pdf')
par(mfrow = c(3,1))
xlim = c(0,1)
hist(pred.unif.congo.bc$mean[best.ix], xlim = xlim)
rug(obs_congo, col = 'red', lwd = 3)
hist(pred.unif.seasia.bc$mean[best.ix], xlim = xlim)
rug(obs_seasia, col = 'red', lwd = 3)
hist(pred.unif.amaz.bc$mean[best.ix], xlim = xlim)
rug(obs_amazon, col = 'red', lwd = 3)
dev.off()

# Which parts of input space are less implausible than the default parameters,
# whenbias corrected to the correct temperature and precipitation?

better.ix.amaz.bc = which(amaz.impl.bc < amaz.impl.bc.default)
better.ix.seasia.bc = which(seasia.impl.bc < seasia.impl.bc.default)
better.ix.congo.bc = which(congo.impl.bc < congo.impl.bc.default)

better.bc.ix = intersect(intersect(better.ix.amaz.bc,better.ix.seasia.bc ), better.ix.congo.bc)

# These are near the edge - might well be uncertainty driving.
pdf(width = 7, height = 7, file = 'graphics/better_bc_default.pdf')
pairs(rbind(X.unif[better.bc.ix, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

# Where do we do better than default parameters?
pred.amaz.bc$mean - obs_amazon

smaller.error.ix.amaz = which(abs(pred.unif.amaz.bc$mean - obs_amazon) < abs(pred.amaz.bc$mean - obs_amazon))
smaller.error.ix.seasia = which(abs(pred.unif.seasia.bc$mean - obs_seasia) < abs(pred.seasia.bc$mean - obs_seasia))
smaller.error.ix.congo = which(abs(pred.unif.congo.bc$mean - obs_congo) < abs(pred.congo.bc$mean - obs_congo))

smaller.bc.ix = intersect(intersect(smaller.error.ix.amaz,smaller.error.ix.seasia ), smaller.error.ix.congo)

(length(smaller.bc.ix) / n) * 100
# These are near the edge - might well be uncertainty driving.
pdf(width = 7, height = 7, file = 'graphics/smaller_error_bc_default.pdf')
pairs(rbind(X.unif[smaller.bc.ix, ], X.stan.norm), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)
dev.off()

pdf(width = 7, height = 7, file = 'graphics/smaller_ae_hists.pdf')
par(mfrow = c(3,1))
xlim = c(0,1)
hist(pred.unif.congo.bc$mean[smaller.error.ix.congo], xlim = xlim)
rug(obs_congo, col = 'red', lwd = 3)
hist(pred.unif.seasia.bc$mean[smaller.error.ix.seasia], xlim = xlim)
rug(obs_seasia, col = 'red', lwd = 3)
hist(pred.unif.amaz.bc$mean[smaller.error.ix.amaz], xlim = xlim)
rug(obs_amazon, col = 'red', lwd = 3)
dev.off()


# --------------------------------------------------------
# What portion of Temperature and Precip space is NROY
# for the default parameters?
# --------------------------------------------------------

# set at default parameters and vary T and P together
# keep points with I < 3
# sample temperature and precip

n = 100000
X.tp = samp.unif(n = n, mins = c(0,0), maxes = c(1,1))
test = matrix()

X.climate = cbind(matrix(rep(X.stan.norm,n), nrow = n, byrow = TRUE), X.tp)
colnames(X.climate) = colnames(X_tropics_norm)

pred.climate = predict(fit.tropics, newdata = X.climate, type = 'UK')

impl.climate.amaz = impl(em = pred.climate$mean, em.sd = pred.climate$sd,
                         disc = disc, obs = obs_amazon,
                         disc.sd = disc.sd,
                         obs.sd = obs.sd)
nroy.ix.climate.amaz = which(impl.climate.amaz < 3)

# South East Asia
impl.climate.seasia = impl(em = pred.climate$mean, em.sd = pred.climate$sd,
                           disc = disc, obs = obs_seasia,
                           disc.sd = disc.sd,
                           obs.sd = obs.sd)
nroy.ix.climate.seasia = which(impl.climate.seasia < 3)

# Central Africa
impl.climate.congo = impl(em = pred.climate$mean, em.sd = pred.climate$sd,
                          disc = disc, obs = obs_congo,
                          disc.sd = disc.sd,
                          obs.sd = obs.sd)
nroy.ix.climate.congo = which(impl.climate.congo < 3)

pdf(width = 7, height = 8, file = 'graphics/nroy_climate.pdf')
par(mfrow = c(2,2), las = 1)

plot(X.climate[nroy.ix.climate.amaz, c(8,9)], type = 'n', xlab = 'Normalised Regional Mean Temperature',
     ylab = 'Normalised Regional Mean  Precipitation',
     main = 'Amazon')
dfunc.up(X.climate[nroy.ix.climate.amaz, 8], X.climate[nroy.ix.climate.amaz, 9])
points(tp.amaz.norm, pch =21, col = 'black', bg = 'red', cex = 1.5)

plot(X.climate[nroy.ix.climate.seasia, c(8,9)], type = 'n', xlab = 'Normalised Regional Mean Temperature',
     ylab = 'Normalised Regional Mean  Precipitation',
     main = 'South East Asia')
dfunc.up(X.climate[nroy.ix.climate.seasia, 8], X.climate[nroy.ix.climate.seasia, 9])
points(tp.seasia.norm, pch =21, col = 'black', bg = 'red', cex = 1.5)

plot(X.climate[nroy.ix.climate.congo, c(8,9)], type = 'n', xlab = 'Normalised Regional Mean Temperature',
     ylab = 'Normalised Regional Mean Precipitation',
     main = 'Central Africa')
dfunc.up(X.climate[nroy.ix.climate.congo, 8], X.climate[nroy.ix.climate.congo, 9])
points(tp.congo.norm, pch =21, col = 'black', bg = 'red', cex = 1.5)
dev.off()


# -------------------------------------------------------------------------------------
# Response to reviewers section:
# What is the uncertainty in the input space that is NROY with respect
# to climates
# NOTE! X.climate only varies across the climate parameters
# -------------------------------------------------------------------------------------

dev.new()
par(mfrow = c(2,2))

cplot(X.climate[nroy.ix.climate.amaz, 8], X.climate[nroy.ix.climate.amaz, 9],
      pred.climate$mean[nroy.ix.climate.amaz],
      cols = blues,
      legend.title =  "emulator mean")

cplot(X.climate[nroy.ix.climate.seasia, 8], X.climate[nroy.ix.climate.seasia, 9],
      pred.climate$mean[nroy.ix.climate.seasia],
      cols = blues,
      legend.title =  "emulator mean")

cplot(X.climate[nroy.ix.climate.congo, 8], X.climate[nroy.ix.climate.congo, 9],
      pred.climate$mean[nroy.ix.climate.congo],
      cols = blues,
      legend.title =  "emulator mean"
      )



dev.new()
par(mfrow = c(2,2))

cplot(X.climate[nroy.ix.climate.amaz, 8], X.climate[nroy.ix.climate.amaz, 9],
      pred.climate$sd[nroy.ix.climate.amaz],
      cols = blues,
      legend.title =  "emulator sd")

cplot(X.climate[nroy.ix.climate.seasia, 8], X.climate[nroy.ix.climate.seasia, 9],
      pred.climate$sd[nroy.ix.climate.seasia],
      cols = blues,
      legend.title =  "emulator sd")

cplot(X.climate[nroy.ix.climate.congo, 8], X.climate[nroy.ix.climate.congo, 9],
      pred.climate$sd[nroy.ix.climate.congo],
      cols = blues,
      legend.title =  "emulator sd")


# ---------------------------------------------------------------------------------
# Response to reviewers
# Monte carlo filtering for sensitivity analysis
# ---------------------------------------------------------------------------------

# Uniform sample from across parameter space
# Split the sample into 'behavioural' (NROY) and 'Non behavioural (Ruled Out)
# Build cdfs of the marginal distributions in each case
# Perform a KS test to see if the smaples are drawn from different distributions
# The KS statistic is an indicator of the importance of the parameter in splitting the
# samples.

# "Not in" function
'%!in%' <- function(x,y)!('%in%'(x,y))

mcf = function(X, nroy.ix){

  ## Monte Carlo Filtering function
  ## X   ............... Complete sample from input space
  ## nroy.ix ........... index of cases of X which are NROY (Not Ruled Out Yet), or 'behavioural'.

  ## produces ks statistic for each column of the input matrix X
  ## A larger ks statistic means that input is more important for
  ## determining if a sample is NROY or not

  X.nroy = X[nroy.ix, ]

  ref = 1:nrow(X)
  ro.ix = which(ref %!in% nroy.ix)
  X.ro = X[ro.ix, ]

  kss = rep(NA, length = ncol(X))
  for(i in 1:ncol(X)){

    ks = ks.test(X.ro[,i], X.nroy[,i])
    kss[i] = ks$statistic

  }

  out = kss
  out
}

# First, use MCF on the design
#
# Calculate the implausibility of each input point for
# each set of observations

run.impl.amaz = impl(em = Y_tropics,
  em.sd = 0,
  disc = 0,
  disc.sd = 0,
  obs = obs_amazon,
  obs.sd = 0.075)

run.nroy.ix.amaz = which(run.impl.amaz < 3)

run.impl.seasia = impl(em = Y_tropics,
  em.sd = 0,
  disc = 0,
  disc.sd = 0,
  obs = obs_seasia,
  obs.sd = 0.075)

run.nroy.ix.seasia = which(run.impl.seasia < 3)

run.impl.congo = impl(em = Y_tropics,
  em.sd = 0,
  disc = 0,
  disc.sd = 0,
  obs = obs_congo,
  obs.sd = 0.075)

run.nroy.ix.congo = which(run.impl.congo < 3)


mcf.amaz = mcf(X_tropics, run.nroy.ix.amaz)
mcf.seasia = mcf(X_tropics, run.nroy.ix.seasia)
mcf.congo = mcf(X_tropics, run.nroy.ix.congo)

dev.new(width = 6, height = 6)
par(mar = c(8,4,3,1))
plot(1:ncol(X_tropics), mcf.amaz, col = col.amaz, pch = 19,
     ylim = c(0,0.4),
     ylab = 'MCF sensitivity', xlab = '',
     axes = FALSE,
     pty = 'n'
     )
abline(v = 1:ncol(X_tropics), lty = 'dashed', col = 'lightgrey')
points(1:ncol(X_tropics), mcf.amaz, col = col.amaz, pch = 19)
points(1:ncol(X_tropics), mcf.seasia, col = col.seasia, pch = 19)
points(1:ncol(X_tropics), mcf.congo, col = col.congo, pch = 19)

axis(side = 1, labels = colnames(X_tropics), las = 2, at = 1:ncol(X_tropics))
axis(side = 2)

# -------------------------------------------------------------------------
# Generate uncertainty estimates on the MCF by emulating and
# bootstrapping the samples.
# -------------------------------------------------------------------------

mcf.emboot = function(X, emfit, bootcol = c(8,9),
  disc, disc.sd, obs, obs.sd, thres = 3, n.mcf = 1000, n.reps = 3000){
  
  ## Function that does Monte Carlo Filtering using an emulated sample.
  ## Inputs for emulation are sampled from the unit cube apart from
  ## those in columns bootcol, which are bootstrapped from the design.
  ##

  em.mcfmat = matrix(nrow = n.reps, ncol = ncol(X))

  for(i in 1:n.reps){

                                        # Sample from uniform distributions for the
                                        # standard input parameters
    X.mcf = samp.unif(n = n.mcf, mins = rep(0, ncol(X)), maxes = rep(1, ncol(X)))
    
                                        # Sample from the model run inputs for the temp and precip
                                        # (here indicated by bootcol columns)
    X.tp.ix = sample(1:nrow(X), size = n.mcf, replace = TRUE)
    X.tp.runsamp = X[X.tp.ix , bootcol]
    
                                        # bind the samples together
    X.mcf[, bootcol] = X.tp.runsamp
    colnames(X.mcf) <- colnames(X)

    # Predict model output at the sampled inputs
    pred.mcf = predict(emfit, newdata = X.mcf, type = 'UK')

    # find the implausibility of the predicted inputs
    em.impl.mcf = impl(em = pred.mcf$mean,
    em.sd = pred.mcf$sd,
      disc = disc,
      disc.sd = disc.sd,
      obs = obs,
      obs.sd = obs.sd)

    # Which part of the sample is NROY (or "behavioural")
    em.nroy.ix = which(em.impl.mcf < thres)
    
    em.mcf= mcf(X.mcf, em.nroy.ix)
    em.mcfmat[i, ] = em.mcf
    
  }
  
  mcf.mean = apply(em.mcfmat, 2, mean)
  mcf.sd = apply(em.mcfmat, 2, sd)


  return(list(mean = mcf.mean, sd = mcf.sd))
}


# How big might the uncertainty bounds be if we use just 300 points
# (as in the ensemble) to estimate the MCF sensitivity analysis indices?

n.mcf.seq = c(seq(from = 100, to = 1000, by = 100), 1500, 2000, 3000)
mcf.seq.mean = matrix(NA, nrow = length(n.mcf.seq), ncol = ncol(X_tropics_norm))
mcf.seq.sd = matrix(NA, nrow = length(n.mcf.seq), ncol = ncol(X_tropics_norm))

for(i in 1:length(n.mcf.seq)){

mcf.em.amaz.seq = mcf.emboot(X = X_tropics_norm, em = tropics_fit,
  bootcol = c(8,9), disc = 0, disc.sd = 0, obs = obs_amazon, obs.sd = 0.05,
  thres = 3, n.mcf = n.mcf.seq[i], n.reps = 1000)

mcf.seq.mean[i, ] = mcf.em.amaz.seq$mean
mcf.seq.sd[i, ] = mcf.em.amaz.seq$sd

}

# What is our estimate of MCF uncertainty when we have 300 ensemble members?
mcf.em.amaz.300 = mcf.emboot(X = X_tropics_norm, em = tropics_fit,
  bootcol = c(8,9), disc = 0, disc.sd = 0, obs = obs_amazon, obs.sd = 0.05,
  thres = 3, n.mcf = 300, n.reps = 1000)

mcf.em.seasia.300 = mcf.emboot(X = X_tropics_norm, em = tropics_fit,
  bootcol = c(8,9), disc = 0, disc.sd = 0, obs = obs_seasia, obs.sd = 0.05,
  thres = 3, n.mcf = 300, n.reps = 1000)

mcf.em.congo.300 = mcf.emboot(X = X_tropics_norm, em = tropics_fit,
  bootcol = c(8,9), disc = 0, disc.sd = 0, obs = obs_congo, obs.sd = 0.05,
  thres = 3, n.mcf = 300, n.reps = 1000)



cbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# The Mean and the Standard deviation of the estimate
# of MCF sensitivity both drop as the number of MCF samples increases.

#dev.new(width = 10, height = 8)
pdf(width = 10, height = 8, file = 'graphics/mcf_mean_sd_vs_n.pdf')
par(mfrow = c(1,2), las = 1)

matplot(n.mcf.seq, mcf.seq.mean, main = 'Mean', xlab = 'Emulated Ensemble members', ylab = 'KS statistic Mean', type = 'o', col = cbPal)

matplot(n.mcf.seq, mcf.seq.sd, main = 'Standard deviation', xlab = 'Emulated Ensemble members', ylab = 'KS statistic standard deviation', type = 'o', col = cbPal)

legend('topright', pch = as.character(1:9), legend = colnames(X_tropics_norm), col = cbPal, text.col = cbPal)
dev.off()


# This puts the mean and estimate MCF sensitivity indices in context with their
# estimated uncertainty.
dev.new(width = 6, height = 10)
matplot(n.mcf.seq, mcf.seq.mean, main = 'Mean', xlab = 'Ensemble members', ylab = 'MCF Sensitivity Index', type = 'o', col = rep(cbPal,2), ylim = c(0,0.35), pch = 19, lwd = 1.2, lty = 'solid')

for(i in 1: ncol(mcf.seq.mean)){
  
arrows(x0 = n.mcf.seq, y0 = mcf.seq.mean[,i] - (mcf.seq.sd[,i]),
         x1 = n.mcf.seq, y1 = mcf.seq.mean[,i] + (mcf.seq.sd[, i]),
         length=0.05, angle=90, code=3, col = rep(cbPal,2)[i],lwd = 1.2
         )
}


# Calculate MCF indices with 5000 emulated ensemble members.
# Bootstrap uncertainty estimates.
mcf.em.amaz = mcf.emboot(X = X_tropics_norm, em = tropics_fit,
  bootcol = c(8,9), disc = 0, disc.sd = 0, obs = obs_amazon, obs.sd = 0.05,
  thres = 3, n.mcf = 5000, n.reps = 1000)

mcf.em.seasia = mcf.emboot(X = X_tropics_norm, em = tropics_fit,
  bootcol = c(8,9), disc = 0, disc.sd = 0, obs = obs_seasia, obs.sd = 0.05,
  thres = 3, n.mcf = 5000, n.reps = 1000)

mcf.em.congo = mcf.emboot(X = X_tropics_norm, em = tropics_fit,
  bootcol = c(8,9), disc = 0, disc.sd = 0, obs = obs_congo, obs.sd = 0.05,
  thres = 3, n.mcf = 5000, n.reps = 1000)



# Plot both the run-generated and emulated MCF sensitivity
#dev.new(width = 8, height = 6)
pdf(file = 'graphics/mcf.pdf', width = 8, height = 6)
par(las = 1, mar = c(8,4,4,2))
ylim = c(0,0.43)
plot((1:length(mcf.em.amaz$mean))-0.2, mcf.em.amaz$mean,
     pch = 19, col = col.amaz, ylim = ylim, xlim = c(0.5,9.5),
     pty = 'n', xaxs = 'i', yaxs = 'i',
     xlab = '', ylab = 'KS statistic',
     axes = FALSE)

i = seq(from = 1, to = 10, by = 2)
rect(i-0.5, ylim[1], i+0.5, ylim[2], col = "lightgrey", border=NA)

points((1:length(mcf.em.amaz$mean))-0.2, mcf.em.amaz$mean, pch = 19, col = col.amaz)

arrows(x0 = (1:length(mcf.em.amaz$mean)) - 0.2, y0 = mcf.em.amaz$mean - (2*mcf.em.amaz$sd ),
         x1 = (1:length(mcf.em.amaz$mean)) - 0.2, y1 = mcf.em.amaz$mean + (2*mcf.em.amaz$sd),
         col = col.amaz, length=0.05, angle=90, code=3)

points((1:length(mcf.em.amaz$mean))-0.2, mcf.amaz, pch = 21, col = col.amaz)

arrows(x0 = 1:length(mcf.em.amaz$mean)-0.2, y0 = mcf.amaz - (2*mcf.em.amaz.300$sd ),
       x1 = 1:length(mcf.em.amaz$mean)-0.2, y1 = mcf.amaz + (2*mcf.em.amaz.300$sd),
       lty = 'dotted',
       col = col.amaz, length=0.05, angle=90, code=3)


points(1:length(mcf.em.seasia$mean), mcf.em.seasia$mean, pch = 19, col = col.seasia)

arrows(x0 = 1:length(mcf.em.seasia$mean), y0 = mcf.em.seasia$mean - (2*mcf.em.seasia$sd ),
         x1 = 1:length(mcf.em.seasia$mean), y1 = mcf.em.seasia$mean + (2*mcf.em.seasia$sd),
         col = col.seasia, length=0.05, angle=90, code=3)


points((1:length(mcf.em.amaz$mean)), mcf.seasia, pch = 21, col = col.seasia)

arrows(x0 = 1:length(mcf.em.amaz$mean), y0 = mcf.seasia - (2*mcf.em.seasia.300$sd ),
       x1 = 1:length(mcf.em.amaz$mean), y1 = mcf.seasia + (2*mcf.em.seasia.300$sd),
       lty = 'dotted',
       col = col.seasia, length=0.05, angle=90, code=3)


points((1:length(mcf.em.congo$mean))+0.2, mcf.em.congo$mean, pch = 19, col = col.congo)

arrows(x0 = (1:length(mcf.em.congo$mean))+0.2, y0 = mcf.em.congo$mean - (2*mcf.em.congo$sd ),
         x1 = (1:length(mcf.em.congo$mean))+0.2, y1 = mcf.em.congo$mean + (2*mcf.em.congo$sd),
         col = col.congo,length=0.05, angle=90, code=3)

points((1:length(mcf.em.congo$mean))+0.2, mcf.congo, pch = 21, col = col.congo)

arrows(x0 = 1:length(mcf.em.amaz$mean)+0.2, y0 = mcf.congo - (2*mcf.em.congo.300$sd ),
       x1 = 1:length(mcf.em.amaz$mean) +0.2, y1 = mcf.congo + (2*mcf.em.congo.300$sd),
       lty = 'dotted',
       col = col.congo, length=0.05, angle=90, code=3)




axis(1, labels = colnames(X_tropics_norm), at = 1:9, las = 2)
axis(2)

legend('topleft',legend = c('Amazon','SE Asia', 'C Africa'),
       col = c(col.amaz, col.seasia, col.congo), pch = 19, bty = 'n')
text(0.5, 0.32, 'Error bars indicate \n \u00B1 2 standard deviations',
     pos  = 4, col = 'black',cex = 0.8 )
text(0.5, 0.29, 'Open points & dotted lines indicate ensemble-only results',
     pos  = 4, col = 'black',cex = 0.8 )
dev.off()


# How does this sensitivity analysis measure up to the FAST99 version?
pdf(width = 12, height = 7, file = 'graphics/fast99_vs_mcf2.pdf')
#dev.new(width = 12, height = 7)
par(mfrow = c(1,2), mar = c(5,5,3,2), las = 1)

plot(print(fast.tell)[,1], mcf.amaz, col = col.amaz, pch = as.character(1:9),
     ylim = c(0,0.42), xlim = c(0,0.32),
     xlab = 'FAST99 first-order sensitivity',
     ylab = 'MCF sensitivity (KS statistic)',
     main = 'MCF using 300 ensemble members',
     pty = 'n'
     )

arrows(x0 = print(fast.tell)[,1], y0 =  mcf.amaz - (2*mcf.em.amaz.300$sd),
         x1 = print(fast.tell)[,1], y1 =  mcf.amaz + (2*mcf.em.amaz.300$sd),
         col = col.amaz, length=0.05, angle=90, code=3,
       lty = 'solid', lwd = 0.8)

arrows(x0 = print(fast.tell)[,1], y0 =  mcf.seasia - (2*mcf.em.seasia.300$sd),
         x1 = print(fast.tell)[,1], y1 =  mcf.seasia + (2*mcf.em.seasia.300$sd),
         col = col.seasia, length=0.05, angle=90, code=3,
       lty = 'solid', lwd = 0.8)

arrows(x0 = print(fast.tell)[,1], y0 =  mcf.congo - (2*mcf.em.congo.300$sd),
         x1 = print(fast.tell)[,1], y1 =  mcf.congo + (2*mcf.em.congo.300$sd),
         col = col.congo, length=0.05, angle=90, code=3,
       lty = 'solid', lwd = 0.8)


points(print(fast.tell)[,1], mcf.amaz, col = col.amaz, pch = as.character(1:9), font = 2)
points(print(fast.tell)[,1], mcf.seasia, col = col.seasia, pch = as.character(1:9), font = 2)
points(print(fast.tell)[,1], mcf.congo, col = col.congo, pch = as.character(1:9), font = 2)

legend('topleft', pch = as.character(1:9), legend = colnames(X_tropics_norm), cex = 0.8, bty = 'n')

legend('top', lty = 'solid', legend = c('Amazon', 'SE Asia', 'C Africa'),
       text.col = c(col.amaz,col.seasia, col.congo),
       col = c(col.amaz,col.seasia, col.congo)
       , cex = 0.8, bty = 'n')

abline(0,1, lty = 'dashed')


plot(print(fast.tell)[,1], mcf.em.amaz$mean, col = col.amaz, pch = as.character(1:9),
     ylim = c(0,0.42), xlim = c(0,0.32),
     xlab = 'FAST99 first-order sensitivity',
     ylab = 'MCF sensitivity (KS statistic)',
     pty = 'n',
     main = 'MCF using 5000 emulated ensemble members'
     )

arrows(x0 = print(fast.tell)[,1], y0 =  mcf.em.amaz$mean - (2*mcf.em.amaz$sd),
         x1 = print(fast.tell)[,1], y1 =  mcf.em.amaz$mean + (2*mcf.em.amaz$sd),
         col = col.amaz,length=0.05, angle=90, code=3)

arrows(x0 = print(fast.tell)[,1], y0 =  mcf.em.seasia$mean - (2*mcf.em.seasia$sd),
         x1 = print(fast.tell)[,1], y1 =  mcf.em.seasia$mean + (2*mcf.em.seasia$sd),
         col = col.seasia,length=0.05, angle=90, code=3)

arrows(x0 = print(fast.tell)[,1], y0 =  mcf.em.congo$mean - (2*mcf.em.congo$sd),
         x1 = print(fast.tell)[,1], y1 =  mcf.em.congo$mean + (2*mcf.em.congo$sd),
         col = col.congo,length=0.05, angle=90, code=3)

points(print(fast.tell)[,1], mcf.em.amaz$mean, col = col.amaz, pch = as.character(1:9), font = 2)
points(print(fast.tell)[,1], mcf.em.seasia$mean, col = col.seasia, pch = as.character(1:9), font = 2)
points(print(fast.tell)[,1], mcf.em.congo$mean, col = col.congo, pch = as.character(1:9), font = 2)


abline(0,1, lty = 'dashed')

dev.off()


# --------------------------------------------------------------
# Analysis suggested by Michael Goldstein - 
# What value does the model add over just using T and P to
# fit the data?
# (not in paper)
# -------------------------------------------------------------

if(run_diagnostics){
  # Produce genuine LOO for all these, put them together and compare with 
  # true.loo
  fit.x.amaz = km(~., design = X.norm, response=famous_agg$AMAZ_MOD_FRAC)
  fit.x.seasia = km(~., design = X.norm, response=famous_agg$SEASIA_MOD_FRAC)
  fit.x.congo = km(~., design = X.norm, response=famous_agg$CONGO_MOD_FRAC)
  
  # This is much quicker than adding them all together!
  true.loo.amaz = true.loo(X = X.norm, y = famous_agg$AMAZ_MOD_FRAC)
  true.loo.seasia = true.loo(X = X.norm, y = famous_agg$SEASIA_MOD_FRAC)
  true.loo.congo = true.loo(X = X.norm, y = famous_agg$CONGO_MOD_FRAC)
  
  true.loo.X.mean = c(true.loo.amaz$mean, true.loo.seasia$mean, true.loo.congo$mean)
  true.loo.X.sd = c(true.loo.amaz$sd, true.loo.seasia$sd, true.loo.congo$sd)
  
  # Mean absolute error is about 0.06 or 6%
  print(paste('Just X mean absolute cross validation error =', mean(abs(true.loo.X.mean - Y_tropics))))
  
  plot(Y_tropics, true.loo.X.mean)
  pdf(width = 6, height = 6, file = 'graphics/true_loo_X.pdf' )
  xlim = c(-0.05, 1.05)
  ylim = c(-0.05, 1.05)
  par(las =1)
  plot(Y_tropics, true.loo.X.mean, pch = 20,
       xlab = 'observation', ylab = 'prediction',
       col = col.tropics,
       xlim = xlim, 
       ylim = ylim,
       bty = 'n',
       axes = FALSE,
       xaxs = 'i', yaxs = 'i')
  
  segments(x0 = Y_tropics, y0 = true.loo.X.mean - (2*true.loo.X.sd),
           x1 = Y_tropics, y1 = true.loo.X.mean +(2*true.loo.X.sd),
           col = col.tropics)
  axis(1, pos = 0, col = 'grey')
  axis(2, pos = 0, col = 'grey')
  abline(0,1, col = 'grey')
  legend('top', legend = c('Amazon', 'Asia', 'Africa'),
         pch = 20, col = c(col.amaz, col.seasia, col.congo),
         bty = 'n')
  dev.off()
  
  
  # Mean absolute error is about 0.03 or 3%
  print(paste('mean absolute cross validation error = ', mean(abs(true.loo.all$mean - Y_tropics))))
  
  # This is much quicker than adding them all together!
  true.loo.tp.amaz = true.loo(X = X_tropics_norm[1:100, 8:9], y = famous_agg$AMAZ_MOD_FRAC)
  true.loo.tp.seasia = true.loo(X = X_tropics_norm[1:100, 8:9], y = famous_agg$SEASIA_MOD_FRAC)
  true.loo.tp.congo = true.loo(X = X_tropics_norm[1:100, 8:9], y = famous_agg$CONGO_MOD_FRAC)
  
  true.loo.tp.mean = c(true.loo.tp.amaz$mean, true.loo.tp.seasia$mean, true.loo.tp.congo$mean)
  true.loo.tp.sd = c(true.loo.tp.amaz$sd, true.loo.tp.seasia$sd, true.loo.tp.congo$sd)
  true.loo.tp.err  = Y_tropics
  
  print(paste('mean absolute cross validation error = ', mean(abs(true.loo.tp.mean - Y_tropics))))
  pdf(width = 6, height = 6, file = 'graphics/true_loo_tp.pdf' )
  xlim = c(-0.05, 1.05)
  ylim = c(-0.05, 1.05)
  par(las =1)
  plot(Y_tropics, true.loo.tp.mean, pch = 20,
       xlab = 'observation', ylab = 'prediction',
       col = col.tropics,
       xlim = xlim, 
       ylim = ylim,
       bty = 'n',
       axes = FALSE,
       xaxs = 'i', yaxs = 'i')
  
  segments(x0 = Y_tropics, y0 = true.loo.X.mean - (2*true.loo.tp.sd),
           x1 = Y_tropics, y1 = true.loo.X.mean +(2*true.loo.tp.sd),
           col = col.tropics)
  axis(1, pos = 0, col = 'grey')
  axis(2, pos = 0, col = 'grey')
  abline(0,1, col = 'grey')
  legend('top', legend = c('Amazon', 'Asia', 'Africa'),
         pch = 20, col = c(col.amaz, col.seasia, col.congo),
         bty = 'n')
  dev.off()
  
  
  # fit using just temperature and precip
  fit.tp  = km(~., design = X_tropics_norm[, 8:9], response=Y_tropics)
  pred.tp = leaveOneOut.km(fit.tp, type="UK", trend.reestim=TRUE)
  
  plot(Y_tropics, pred.tp$mean )
  fit.tp.resid = pred.tp$mean - Y_tropics
  
  # Have to split these out per-forest
  fit.resid.amazon = km(~., design = X.norm, response=fit.tp.resid[1:100])
  fit.resid.seasia = km(~., design = X.norm, response=fit.tp.resid[101:200])
  fit.resid.congo = km(~., design = X.norm, response=fit.tp.resid[201:300])
  
  plot(fit.resid.congo)
  
}else{print('skipping diagnostics')}

# ---------------------------------------------------------------------------
# Plot maps of the tropical forests.
# ---------------------------------------------------------------------------
remap.famous = function(dat,longs,lats, shift = FALSE){
  # reshape a map in vector form so that fields() package function image.plot() 
  #  (for example) will plot it correctly
  mat = matrix(dat, nrow=length(longs), ncol=length(lats))[ ,length(lats):1]
  if(shift){
    block1.ix = which(longs <= shift)
    block2.ix = which(longs > shift)
    mat.shift = rbind(mat[ block2.ix, ], mat[block1.ix, ]) 
    out = mat.shift
  }
  else{
    out = mat
  }
  out
}

remap.tropics = function(dat,lats, upper, lower){
  ix = which(lats <= upper & lats >=lower)
  out = dat[,ix]
}

blockswap = function(dat, longs,lats, shift){
  # takes a map like matrix already in the right format
  # for mapping in image.plot, and plots it from a different
  # longitude
  block1.ix = which(longs <= shift)
  block2.ix = which(longs > shift)
  mat.shift = rbind(dat[ block2.ix, ], dat[block1.ix, ]) 
  out = mat.shift
}


famous.example = blockswap(remap.famous(bl.frac.ens[1,], longs = longs, lats = lats),
                           longs = longs, lats = lats, shift = 180)


load('forest_fraction_obs_map.RData')
# HadGEM2 family resolution
#lats = 1.25 * 1.875
#bl.obs.dat = read.table('forest_fraction_obs_map_v2.txt', na.strings = '-1.073741824000000000e+09')

obslats = seq(from = -90, to = 90, length.out =  dim(bl.obs.dat)[2])
obslongs = seq(from = 0, to = (360-1.875), by = 1.875)

bl.obs.map = blockswap(t(as.matrix(bl.obs.dat)), longs = obslongs, lats = obslats, shift = 180)

#bl.dat.regrid = read.table('forest_fraction_obs_map_regrid_v2.txt', na.strings = '-1.073741824000000000e+09')
bl.obs.map.regrid = blockswap(t(as.matrix(bl.dat.regrid)), longs = longs, lats = lats, shift = 180)

pdf(width = 5, height = 8, file = 'graphics/map_comparison.pdf' )
par(bg = 'lightgrey', mfrow = c(2,1), oma = c(4,0,0,0), mar = c(4,1,3,1))
image(bl.obs.map, col = yg, zlim = c(0,1),  axes = FALSE, main = 'Observations')
image(bl.obs.map.regrid, col = yg, zlim = c(0,1),  axes = FALSE, main = 'Regridded observations')
reset()
par(oma = c(1,0,0,0))
image.plot(bl.obs.map, zlim = c(0,1), legend.only = TRUE, horizontal = TRUE, 
           col = yg, legend.shrink = 0.6, legend.width = 0.7,
           legend.lab = 'Broadleaf forest fraction')
dev.off()

blmeans = rep(NA, 100)
for(i in 1:100){
  blmeans[i] = mean(bl.frac.ens[i,], na.rm = TRUE)
}

bl.ix = order(blmeans)

pdf(file = 'graphics/tropics_maps_yg.pdf', width = 8, height = 8)
par(mfrow = c(13,8), mar = c(0.2, 0.2, 0.2, 0.2), bg = 'lightgrey',
    oma = c(9,0.2,0.2,0.2))

for(i in bl.ix){
  
  map = remap.famous(bl.frac.ens[i,], longs = longs, lats = lats, shift = 180)
  test.trop = remap.tropics(map,lats = rev(lats), upper = 60, lower = -60)
  image(test.trop, axes = FALSE, col = yg, zlim = c(0,1))
  
}

reset()
par(oma = c(1,0,0,0))
image.plot(test.trop, zlim = c(0,1), legend.only = TRUE, horizontal = TRUE, 
           col = yg, legend.shrink = 0.6, legend.width = 0.7,
           legend.lab = 'Broadleaf Forest fraction')
dev.off()


pdf(file = 'graphics/tropics_anom_maps.pdf', width = 8, height = 8)
par(mfrow = c(13,8), mar = c(0.2, 0.2, 0.2, 0.2), bg = 'lightgrey',
    oma = c(9,0.2,0.2,0.2))

for(i in bl.ix){
  
  bl = remap.famous(bl.frac.ens[i,], longs = longs, lats = lats, shift = 180)
  anom = bl - bl.obs.map.regrid
  image(remap.tropics(anom, lats = rev(lats), upper = 60, lower = -60), axes = FALSE, col = rev(byr), zlim = c(-1,1))
  
}

reset()
par(oma = c(1,0,0,0))
image.plot(anom, zlim = c(-1,1), legend.only = TRUE, horizontal = TRUE, 
           col = rev(byr), legend.shrink = 0.6, legend.width = 0.7,
           legend.lab = 'Broadleaf forest fraction anomaly')
dev.off()




# --------------------------------------------------------------------------------------
# Validation section for response to reviewers
# --------------------------------------------------------------------------------------

# Find the ensemble members nearest the observed temperature and precipitation
# for the three forests, and those at the edges of the temperature and precipitation
# data set. Hold out those members and check how well the GP does in prediction.

# Choose ensemble members near to the edges of the
# temperature and precip parameter space

holdout1.ix <- which(X_tropics_norm[,'MOD_TEMP'] < 0.4 & X_tropics_norm[,'MOD_PRECIP'] < 0.2)

X.test.holdout1 = X_tropics_norm[holdout1.ix, ]
X.train.holdout1 = X_tropics_norm[-holdout1.ix, ]

y.test.holdout1 <- Y_tropics[holdout1.ix]
y.train.holdout1 <- Y_tropics[-holdout1.ix]


pch.tropics = c(rep(pch.amaz, 100), rep(pch.seasia, 100), rep(pch.congo, 100))

#dev.new()
pdf(file = 'graphics/holdout1_location.pdf')

plot(X.test.holdout1[,'MOD_TEMP'], X.test.holdout1[,'MOD_PRECIP'], xlim = c(0,1),
     ylim = c(0,1), pch = pch.tropics[holdout1.ix], col = 'red',
     xlab = 'normalised temperature', ylab = 'normalised preipitation')
points(X.train.holdout1[,'MOD_TEMP'],X.train.holdout1[,'MOD_PRECIP'] , pch =  pch.tropics[-holdout1.ix], col = 'black')
points(tp.amaz.norm, col = 'black', pch = pch.amaz, cex = 2, bg = col.amaz, lwd = 2.5)
points(tp.congo.norm, col = 'black', pch = pch.congo, cex = 2, bg = col.congo, lwd = 2.5)
points(tp.seasia.norm, col = 'black', pch = pch.seasia, cex = 2, bg = col.seasia, lwd = 2.5)

text(tp.amaz.norm, 'Amazon', pos = 4, font = 2)
text(tp.congo.norm, 'Central Africa', pos = 4, font = 2)
text(tp.seasia.norm, 'SE Asia', pos = 4, font = 2)
dev.off()



fit.holdout1 = km(~., design = X.train.holdout1, response = y.train.holdout1)
pred.holdout1 = predict(fit.holdout1, newdata = X.test.holdout1, type = 'UK')

plot(1:6, y.test.holdout1, ylim = c(0,1), pch = 19)
points(1:6, pred.holdout1$mean, ylim = c(0,1), col = 'darkgrey')
segments(x0 = 1:6, y0 = pred.holdout1$mean - (2*pred.holdout1$sd),
         x1 = 1:6, y1 = pred.holdout1$mean + (2*pred.holdout1$sd),
         col = 'darkgrey'
         )


# Sort predictions by forest fraction magnitude for plotting
frac.sort = sort(Y_tropics, index.return = TRUE)

pdf(file = 'graphics/holdout1_vs_loo.pdf', width = 16, height = 4)
# Choose ensemble members near to the edges of the
# temperature and precip parameter space
  
plot(1:300, Y_tropics[frac.sort$ix], pch = 20, xlim = c(-1, 310),ylim = c(-0.1,1), xaxs = 'i',
     xlab = 'ensemble member rank',
     ylab = 'forest fraction')
points(1:300, true.loo.all$mean[frac.sort$ix], pch = 20, col = 'darkgrey')

y0 = (true.loo.all$mean - (2*true.loo.all$sd))[frac.sort$ix]
y1 = (true.loo.all$mean + (2*true.loo.all$sd))[frac.sort$ix]

segments(x0 = 1:300, y0 = y0,
         x1 = 1:300,y1 = y1,
         col = 'darkgrey'
         )

points(302:307, y.test.holdout1, pch = 20, col = 'red')
points(302:307, pred.holdout1$mean, ylim = c(0,1),pch = 20, col = 'darkgrey')

segments(x0 = 302:307, y0 = pred.holdout1$mean - (2*pred.holdout1$sd),
         x1 = 302:307, y1 = pred.holdout1$mean + (2*pred.holdout1$sd),
         col = 'darkgrey'
         )
abline(v = 301, lty = 'dashed')

legend('topleft', pch = 20, col = c('black', 'red', 'darkgrey'),
       legend = c('Ensemble member', 'hold-out sample', 'prediction'), bty = 'n', cex = 0.8)

dev.off()


# Histogram of leave-one-out and holdout errors
pdf(file = 'graphics/holdout1_error_hist.pdf', width = 6, height = 5)
hist(true.loo.all$mean - Y_tropics, col = 'darkgrey',
     xlab = 'prediction error',
     main = '')

legend('topleft', pch = c('|',NA), col = c('red', NA),
       legend = c('hold-out sample','leave-one-out') , cex = 0.8, fill  =c(NA, 'darkgrey'),
       border = c(NA, 'black'))

rug(pred.holdout1$mean - y.test.holdout1, col = 'red', lwd = 2)
dev.off()


# Compare leave-one-out error against hold-out error for the 6 members.


#dev.new(width= 7, height = 4)
pdf(file = 'graphics/loo_v_holdout1_prediction_error.pdf', width = 7, height = 4)
par( las = 1)
plot(1:6, Y_tropics[holdout1.ix], ylim = c(0,1), pch = 19, xlim = c(1, 6.2),
     ylab = 'forest fraction',
     xlab = 'held-out ensemble member')

abline(h = range(Y_tropics), col = 'darkgrey', lty = 'dashed') 

points(1.1:6.1, true.loo.all$mean[holdout1.ix], pch = 19, col = 'blue')

segments(x0 = 1.1:6.1, y0 =  (true.loo.all$mean[holdout1.ix] - 2*true.loo.all$sd[holdout1.ix]),
         x1 = 1.1:6.1, y1 = (true.loo.all$mean[holdout1.ix] + 2*true.loo.all$sd[holdout1.ix]),
         col = 'blue'
         )


points(1.2:6.2, pred.holdout1$mean, col = 'red', pch = 19)
segments(x0 = 1.2:6.2, y0 =  (pred.holdout1$mean - 2*pred.holdout1$sd),
         x1 = 1.2:6.2, y1 = (pred.holdout1$mean + 2*pred.holdout1$sd),
         col = 'red'
         )
legend('bottomright',
       col = c('black', 'blue', 'red'),
       pch = 19,
       legend = c('Ensemble member','LOO prediction', 'holdout prediction'),
       bty = 'n',
       cex = 0.8, pt.cex = 1,
       inset = 0.08
       )
text(1,0.1, labels = 'Dashed lines are ensemble limits', col = 'darkgrey', cex = 0.8,
     pos = 4)
graphics.off()


loo.err.holdout1 = (true.loo.all$mean - Y_tropics)[holdout1.ix]
err.holdout1 = pred.holdout1$mean - y.test.holdout1
dev.new()
plot(1:6, abs(loo.err.holdout1), pch = 19, col = 'blue', ylim = c(0,0.06))
points(1:6, abs(err.holdout1), pch = 19, col = 'red')


# -------------------------------------------------------------------------
# Reviewer 1 is concerned about the impact of a lack of samples in  
# parts of temperature and precipitation space, and how that might impact
# on sensitivity analysis and history matching conclusions
# -------------------------------------------------------------------------


# Code from here is just copied from elsewhere at the moment, needs adjusting
# to apply to these data.

constrained.oaat = function(X, Y, n.oaat = 21, mins, maxes, hold = NULL,
                            predtype = 'UK',
                            nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL,
                            seed=NULL, trace=FALSE, maxit=100,
                            REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, 
                            popsize=100)
  {
  # X        ...   Experiment design
  # Y        ...   Matrix of model output, with variables in columns
  # maxes    ...   Vector maximum tolerable value of variables corresponding
  #                to columns of Y
  # mins     ...   Vector minimum tolerable value of variables corresponding
  #                to columns of Y
  
  # generate oaat design
  d = ncol(X)
  X.norm = normalize(X)
  X.oaat = oaat.design(X.norm, n = n.oaat, med = TRUE, hold = hold)
  colnames(X.oaat) = colnames(X)
  
  # generate ncol(Y) emulators
  p = ncol(Y)
  pred.mean = matrix(NA, ncol = p, nrow = nrow(X.oaat))
  pred.sd = matrix(NA, ncol = p, nrow = nrow(X.oaat))
  
  for(i in 1:p){
    
    y = Y[, i]
    
    em = twoStep.glmnet(X=X.norm, y=y, nugget=nugget, nuggetEstim=nuggetEstim,
                                  noiseVar=noiseVar,
                                  seed=seed, trace=trace, maxit=maxit,
                                 REPORT=REPORT, factr=factr, pgtol=pgtol,
                                 parinit=parinit, popsize=popsize)
    
    oaat.pred = predict(em$emulator, newdata = X.oaat, type = predtype)
    
    # produce the whole oaat emulator output
    pred.mean[, i] = oaat.pred$mean
    pred.sd[, i] = oaat.pred$sd
    
  }
  
  ix.kept = apply(pred.mean, 1, FUN = allin, mins = mins, maxes = maxes)
  
  # Replace out-of-bound rows in X with NA, to make plotting
  # the final results easier
  X.oaat.constr = X.oaat
  X.oaat.constr[ix.kept==FALSE, ] <- NA
  
  pred.constr = pred.mean
  pred.constr[ix.kept==FALSE, ] <- NA
  
  pred.sd.constr = pred.sd
  pred.sd.constr[ix.kept==FALSE, ] <- NA
  
  # keep only the oaat emulator output within constraints
  return(list(X.oaat.constr = X.oaat.constr,
              ix.kept = ix.kept,
              pred.constr = pred.constr,
              pred.sd.constr = pred.sd.constr,
              X.oaat = X.oaat,
              pred.mean = pred.mean,
              pred.sd = pred.sd)
         )
}

inputs.set <- function(X, y, thres, obs, obs.sd = 0, disc = 0, disc.sd = 0, n = 100000, abt = FALSE){ 
  # find a set of inputs that are consistent with a particular
  # set of implausibility (either below or above)
  
  X.mins <- apply(X,2,min)
  X.maxes <- apply(X,2,max)
  X.unif <- samp.unif(n, mins = X.mins, maxes = X.maxes)
  colnames(X.unif) <- colnames(X)
  
  fit <- km(~., design = X, response = y, control = list(trace = FALSE))
  pred <- predict(fit, newdata = X.unif, type = 'UK')
  pred.impl <- impl(em = pred$mean, em.sd = pred$sd,
                    disc = disc, obs = obs, disc.sd = disc.sd, obs.sd = obs.sd)
  
  if(abt){
    # choose those above the threshold 
    ix.bt <- pred.impl > thres
  }
  
  else{
    ix.bt <- pred.impl < thres
  }
  
  X.out <- X.unif[ix.bt, ]
  
  return(list(X.out = X.out, fit = fit, X.unif = X.unif, pred = pred,pred.impl = pred.impl))   
}


plausible.amazon.bc <- inputs.set(X = X_tropics_norm, y = Y_tropics,thres = 3,
                                  obs = obs_amazon,
                                  obs.sd = 0,
                                  disc = 0,
                                  disc.sd = 0.01,
                                  n = 100000,
                                  abt = FALSE)



# Combine the above two to create a sensitivity analysis that excludes parts of input
# parameter space that are deemed ruled out

historymatch.oaat <- function(X, Y,
                              thres,
                              obs,
                              obs.sd = 0, disc = 0, disc.sd = 0,
                              n.oaat = 21, abt = FALSE,
                              hold = NULL,
                              predtype = 'UK',
                              nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL,
                              seed=NULL, trace=FALSE, maxit=100,
                              REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, 
                              popsize=10
                              ){
 # X   ...   design
 # Y   ...   Matrix of outputs, with rows matching rows of X
  


  d = ncol(X)
  X.norm = normalize(X)
  X.oaat = oaat.design(X.norm, n = n.oaat, med = TRUE, hold = hold)
  colnames(X.oaat) = colnames(X)
  
  # generate ncol(Y) emulators
  p = ncol(Y)
  pred.mean = matrix(NA, ncol = p, nrow = nrow(X.oaat))
  pred.sd = matrix(NA, ncol = p, nrow = nrow(X.oaat))
  pred.impl = matrix(NA, ncol = p, nrow = nrow(X.oaat))
  
  for(i in 1:p){
    
    y = Y[, i]
    
    em = twoStep.glmnet(X=X.norm, y=y, nugget=nugget, nuggetEstim=nuggetEstim,
                                  noiseVar=noiseVar,
                                  seed=seed, trace=trace, maxit=maxit,
                                 REPORT=REPORT, factr=factr, pgtol=pgtol,
                                 parinit=parinit, popsize=popsize)
    
    oaat.pred = predict(em$emulator, newdata = X.oaat, type = predtype)
    
    # produce the whole oaat emulator output
    pred.mean[, i] = oaat.pred$mean
    pred.sd[, i] = oaat.pred$sd

    # fix this up for each Y
    pred.impl[, i] = impl(em = pred$mean, em.sd = pred$sd,
                    disc = disc, obs = obs, disc.sd = disc.sd, obs.sd = obs.sd)

  }
  
  # Replace out-of-bound rows in X with NA, to make plotting
  # the final results easier


  max.impl <- apply(pred.impl,2,max)
  ix.bt <- max.impl < thres

  X.oaat.constr = X.oaat
  X.oaat.constr[ix.bt==FALSE, ] <- NA
  
  return(list(X.out = X.out,
              X.oaat = X.oaat,
              pred.impl = pred.impl))  
 
}

# Its going to be tricky - need to think about how the data is put together.
# In the augmented emulator, everything is combined (i.e. there are three different
# observations. Need the emulator from the combined data, but History matching
# individual data streams.

  



normalize.na = function(X, wrt = NULL){ 
  
  f <- function(X){
    (X-min(X, na.rm = TRUE))/(max(X, na.rm = TRUE)-min(X, na.rm = TRUE))
  }
  
  # test to see if we have a matrix, array or data frame
  if(length(dim(X))==2){
    out <- apply(X,2,f)
  }
  
  else{	
    out <- f(X)
  }
  
  if(is.null(wrt) == FALSE){
    # if argument wrt is given
    
    n <- nrow(X)
    mmins <- t(kronecker(apply(wrt,2,min, na.rm = TRUE),t(rep(1,n))))
    mmaxs <- t(kronecker(apply(wrt,2,max, na.rm = TRUE),t(rep(1,n))))
    
    out <- (X-mmins)/(mmaxs-mmins)
    
  }
  
  out
}

# Express the standard in terms of lhs, then normalised matrices
X.stan = c(rep(1, d-1),0)
lhs.range = round( apply(lhs,2,range),1)

# Standard parameters when compared to the initial latin hypercube
X.stan.wrt.lhs = normalize(matrix(X.stan, nrow = 1), wrt = lhs.range)

# Normalize BEFORE putting it in to the SA
# Keep everything in relation to original design


# Need to normalize the constraints too
# Normalize everything compared to the initial data (dat.norm)
mins.norm = normalize.na(matrix(mins.constr, nrow = 1), wrt = dat.norm)
maxes.norm = normalize.na(matrix(maxes.constr, nrow = 1), wrt = dat.norm)


Y.norm = normalize.na(dat.level0, wrt = dat.norm)

glob.const.oaat = constrained.oaat(X = X.level0,
  Y = Y.norm,
  n.oaat = 21,
  mins = mins.norm,
  maxes = maxes.norm,
  hold = X.stan.wrt.lhs
  )


# It's not constraining at the moment: why not?
# Its the corners that are ruled out!!


n = 21
# Here is the full sensitivity analysis

linecols.ext = c('black', paired)
ylim = c(0,1)
pdf(file = 'graphics/global_not_constrained_oaat.pdf', width = 9, height = 9)
#dev.new(width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))
ndat = ncol(glob.const.oaat$pred.mean)
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  
  y.oaat = glob.const.oaat$pred.mean

  plot(glob.const.oaat$X.oaat[ix,i], y.oaat[ix],
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:ndat){
    y.oaat = glob.const.oaat$pred.mean[ix,j]
    lines(glob.const.oaat$X.oaat[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  
}
reset()
legend('top',
       legend = colnames(dat.norm), 
       col = linecols.ext,
       cex = 0.8,
       lwd = 2,
       horiz = TRUE)
dev.off()


linecols.ext = c('black', paired)
ylim = c(0,1)
pdf(file = 'graphics/global_constrained_oaat.pdf', width = 9, height = 9)
#dev.new(width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))
ndat = ncol(glob.const.oaat$pred.constr)
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  
  y.oaat = glob.const.oaat$pred.constr
  
  plot(c(0,1), c(0,1),
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:ndat){
    y.oaat = glob.const.oaat$pred.constr[ix,j]
    lines(glob.const.oaat$X.oaat.constr[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  
}
reset()
legend('top',
       legend = colnames(dat.level0), 
       col = linecols.ext,
       cex = 0.8,
       lwd = 2,
       horiz = TRUE)

dev.off()

