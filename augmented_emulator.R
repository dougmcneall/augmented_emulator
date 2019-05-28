# augmented_emulator.R
# This code augments a Gaussian process emulator using temperature
# and precipitation (usually considered model outputs) as inputs.
# It allows biases in those fields to be accounted for when choosing
# plausible values of input parameters for the land surface component
# of the reduced-resolution climate model, FAMOUS.
# This is the main analysis code for:
# McNeall, D.J., Williams J., Betts, R.A. , Booth, B.B.B., Challenor, P.G. , Good, P.
# & Wiltshire A. (2019) submitted to Geoscientific Model Development
# Contact Doug McNeall dougmcneall@gmail.com @dougmcneall

# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------
library(DiceKriging)
library(RColorBrewer)
library(MASS)
library(fields)
library(parallel)
library(viridisLite)

load('famous_forest_fraction.RData')
load('famous_agg.RData')

# Load specific versions of a github repository
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/emtools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/imptools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/vistools.R')

#source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R')
#source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R')
#source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R')

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

dfunc.up <- function(x,y,...){
  require(MASS)
  require(RColorBrewer)
  
  br <- brewer.pal(9, 'Blues')
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde <- kde2d(x,y)
  image(kde, col = br, add = TRUE)
}

dfunc.up.truth = function(x,y, ...){
  require(MASS)
  require(RColorBrewer)
  
  xtrue <- tail(x,1)
  ytrue <- tail(y,1)
  
  xdash <- head(x, -1)
  ydash <- head(y, -1)
  
  br <- brewer.pal(9, 'Blues')
  # function for plotting 2d kernel density estimates in pairs() plot.
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

# Plot colour as the third dimension, to compare model runs with
# mean emulated surface.
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
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}


X = famous_agg[, 2:8]
X.norm = normalize(X)
X.stan.norm <- normalize(matrix(X.standard, nrow = 1), wrt = X)
colnames(X.stan.norm) = colnames(X)

# Matrix of (repeated) input parameters, modelled temperatures,
# precipitation and dorest fractions.

X_AMAZON = famous_agg[ ,c(2,3,4,5,6,7,8,15,19)] # AMAZON
X_SEASIA = famous_agg[, c(2,3,4,5,6,7,8, 17,21)] # SEASIA
X_CONGO = famous_agg[, c(2,3,4,5,6,7,8,16,20)] # CONGO

colnames(X_AMAZON) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')
colnames(X_SEASIA) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')
colnames(X_CONGO) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')

X_tropics = rbind(X_AMAZON, X_SEASIA, X_CONGO)
X_tropics_norm = normalize(X_tropics)

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


#fit.namerica <- km(~., design=X.norm, response=NAMERICA_MOD_FRAC)
#fit.global <- km(~., design=X.norm, response=GLOB_MOD_FRAC)

standard.amazon <- predict(fit.amazon, newdata=X.stan.norm, type='UK')
standard.congo <- predict(fit.congo, newdata=X.stan.norm, type='UK')
standard.seasia <- predict(fit.seasia, newdata=X.stan.norm, type='UK')

#standard.namerica <- predict(fit.namerica, newdata=X.stan.norm, type='UK')
#standard.global <- predict(fit.global, newdata=X.stan.norm, type='UK')

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

#plot(X_tropics[, 8], X_tropics[, 9], col = 'black', bg = zcolor, pch = 21, cex = 2)

# ----------------------------------------------------------------------
# This plot Deprecated in favour of a quilt plot (below) for the paper 
# ----------------------------------------------------------------------

# Normalize the colours to the background
# the first 300 points are Y_tropics
# ncols = 7
# allz = c(Y_tropics,obs_amazon,obs_seasia, obs_congo, plausible.amazon.bc$pred$mean)
# zcolor = col3rd(n=ncols, pal=viridis(ncols), z = allz) 
# 
# pdf(file = 'graphics/emulated_fraction_vs_temp_precip_pcolcor.pdf',width = 7, height = 7)
# par(las = 1)
# quilt.plot(plausible.amazon.bc$X.unif[,8],
#            plausible.amazon.bc$X.unif[, 9], plausible.amazon.bc$pred$mean,
#            col = viridis(ncols),
#            xlab = 'Normalised regional temperature', 
#            ylab = 'Normalised regional precipitation', 
#            legend.args = list(text = "forest\nfraction",
#                               col="black", cex=1.2, side=3, line=1))
# cex = 1.4
# lwd = 1.5
# points(X_tropics_norm[1:100,8], X_tropics_norm[1:100,9], 
#        col = 'black', bg = zcolor[1:100], pch = 21, cex = cex, lwd = lwd)
# points(X_tropics_norm[101:200,8], X_tropics_norm[101:200,9], 
#        col = 'black', bg = zcolor[101:200], pch = 22, cex = cex, lwd = lwd)
# points(X_tropics_norm[201:300,8], X_tropics_norm[201:300,9], 
#        col = 'black', bg = zcolor[201:300], pch = 24, cex = cex, lwd = lwd)
# 
# points(tp.amaz.norm, col = 'black', pch = 21, cex = 2.5, bg = zcolor[301], lwd = 2)
# points(tp.seasia.norm, col = 'black', pch = 22, cex = 2.5, bg = zcolor[302], lwd = 2)
# points(tp.congo.norm, col = 'black', pch = 24, cex = 2.5, bg = zcolor[303], lwd = 2)
# 
# shadowtext(tp.amaz.norm[1],tp.amaz.norm[2], 'Amazon', pos = 4, font = 2,r =0.2)
# shadowtext(tp.seasia.norm[1], tp.seasia.norm[2], 'SE Asia', pos = 4, font = 2, r = 0.2)
# shadowtext(tp.congo.norm[1], tp.congo.norm[2], 'Central Africa', pos = 4, font = 2, r = 0.2)
# dev.off()

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


# --------------------------------------------------------------
# Analysis suggested by Michael Goldstein - 
# What value does the model add over just using T and P to
# fit the data?
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





