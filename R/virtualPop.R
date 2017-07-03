#' @title Virtual fish population
#'
#' @param K.mu mean K (growth parameter from von Bertalanffy growth function)
#' @param K.cv coefficient of variation on K
#' @param Linf.mu mean Linf (infinite length parameter from von Bertalanffy growth function)
#' @param Linf.cv coefficient of variation on Linf
#' @param ts summer point (range 0 to 1) (parameter from seasonally oscillating von Bertalanffy growth function)
#' @param C strength of seasonal oscillation (range 0 to 1) (parameter from seasonally oscillating von Bertalanffy growth function)
#' @param LWa length-weight relationship constant 'a' (W = a*L^b). Model assumed length in cm and weight in kg.
#' @param LWb length-weight relationship constant 'b' (W = a*L^b). Model assumed length in cm and weight in kg.
#' @param Lmat length at maturity (where 50\% of individuals are mature)
#' @param wmat width between 25\% and 75\% quantiles for Lmat
#' @param rmaxBH parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param betaBH parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param srr.cv coefficient of variation on number of recruits
#' @param repro_wt weight of reproduction (vector of monthly reproduction weight)
#' @param M natural mortality
#' @param Etf  effort (E = F / q); single numeric, numeric vector for effort per year, or matrix for different fleets (columns) and different years (rows)
#' @param qtf catchability (default 0.005); single numeric, numeric vector for effort per year, or matrix for different fleets (columns)  and different years (rows)
#' @param harvest_rate Fishing mortality (i.e. 'F' = C/B); if NaN Etf and qtf are used to estimate the harvest_rate
#' @param L50 minimum length of capture (in cm). Where selectivity equals 0.5. Assumes logistic ogive typical of trawl net selectivity.
#' @param wqs width of selectivity ogive (in cm)
#' @param bin.size resulting bin size for length frequencies (in cm)
#' @param timemin time at start of simulation (in years). Typically set to zero.
#' @param timemax time at end of simulation (in years).
#' @param timemin.date date corresponding to timemin (of "Date" class)
#' @param tincr time increment for simulation (default = 1/12; i.e. 1 month)
#' @param N0 starting number of individuals
#' @param fished_t times when stock is fished; when NA no exploitation simulated
#' @param lfqFrac fraction of fished stock that are sampled for length frequency data (default = 0.1).
#' @param progressBar Logical. Should progress bar be shown in console (Default=TRUE)
#'
#' @description See \code{\link[fishdynr]{dt_growth_soVB}} for information on growth function.
#' The model creates variation in growth based on a mean phi prime value for the population,
#' which describes relationship between individual Linf and K values. See Vakily (1992)
#' for more details.
#'
#' @details The model takes around 5 to 10 years to reach equilibrium, i.e. no biomass changes independent from fishing activity, the actual time is dependent on N0, K.mu, Lmat, repro_wt  and rmax.BH. For the estimation of carrying capacity the first 10 years of the simulation are disregarded and only subsequent years where no fishing took place are used to estimate the annual mean carrying capacity (K). If fishing is simulated for all years or fishing activities start before ten years after simulation start no carrying capacity is estimated.
#'
#' @return a list containing growth parameters and length frequency object
#'
#' @references
#' Vakily, J.M., 1992. Determination and comparison of bivalve growth,
#' with emphasis on Thailand and other tropical areas. WorldFish.
#'
#' Munro, J.L., Pauly, D., 1983. A simple method for comparing the growth
#' of fishes and invertebrates. Fishbyte 1, 5-6.
#'
#' Pauly, D., Munro, J., 1984. Once more on the comparison of growth
#' in fish and invertebrates. Fishbyte (Philippines).
#'
#' @importFrom graphics hist
#' @importFrom stats rlnorm runif weighted.mean
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats qnorm rnorm
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' res <- virtualPop()
#' names(res)
#'
#' op <- par(mfcol=c(2,1), mar=c(4,4,1,1))
#' plot(N ~ dates, data=res$pop, t="l")
#' plot(B ~ dates, data=res$pop, t="l", ylab="B, SSB")
#' lines(SSB ~ dates, data=res$pop, t="l", lty=2)
#' par(op)
#'
#' pal <- colorRampPalette(c("grey30",5,7,2), bias=2)
#' with(res$lfqbin, image(x=dates, y=midLengths, z=t(catch), col=pal(100)))
#'
#' ### biased results with single monthly sample
#' inds <- res$inds[[1]]
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' fit <- glm(mat ~ L, data = inds, family = binomial(link = "logit"))
#' summary(fit)
#'
#' newdat <- data.frame(L = seq(min(inds$L), max(inds$L), length.out=100))
#' newdat$pmat <- pmat_w(newdat$L, Lmat = 40, wmat=40*0.2)
#' pred <- predict(fit, newdata=newdat, se.fit=TRUE)
#' # Combine the hypothetical data and predicted values
#' newdat <- cbind(newdat, pred)
#' # Calculate confidence intervals
#' std <- qnorm(0.95 / 2 + 0.5)
#' newdat$ymin <- fit$family$linkinv(newdat$fit - std * newdat$se.fit)
#' newdat$ymax <- fit$family$linkinv(newdat$fit + std * newdat$se.fit)
#' newdat$fit <- fit$family$linkinv(newdat$fit)  # Rescale to 0-1
#'
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' lines(pmat ~ L, newdat, col=8, lty=3)
#' polygon(
#'   x = c(newdat$L, rev(newdat$L)),
#'   y = c(newdat$ymax, rev(newdat$ymin)),
#'   col = adjustcolor(2, alpha.f = 0.3),
#'   border = adjustcolor(2, alpha.f = 0.3)
#' )
#' lines(fit ~ L, newdat, col=2)
#'
#' lrPerc <- function(alpha, beta, p) (log(p/(1-p))-alpha)/beta
#' ( L50 <- lrPerc(alpha=coef(fit)[1], beta=coef(fit)[2], p=0.5) )
#' lines(x=c(L50,L50,0), y=c(-100,0.5,0.5), lty=2, col=2)
#' text(x=L50, y=0.5, labels = paste0("L50 = ", round(L50,2)), pos=4, col=2 )
#'
#'
#'
#' ### all samples combined
#' inds <- do.call("rbind", res$inds)
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' fit <- glm(mat ~ L, data = inds, family = binomial(link = "logit"))
#' summary(fit)
#'
#' newdat <- data.frame(L = seq(min(inds$L), max(inds$L), length.out=100))
#' newdat$pmat <- pmat_w(newdat$L, Lmat = 40, wmat=40*0.2)
#' pred <- predict(fit, newdata=newdat, se.fit=TRUE)
#' # Combine the hypothetical data and predicted values
#' newdat <- cbind(newdat, pred)
#' # Calculate confidence intervals
#' std <- qnorm(0.95 / 2 + 0.5)
#' newdat$ymin <- fit$family$linkinv(newdat$fit - std * newdat$se.fit)
#' newdat$ymax <- fit$family$linkinv(newdat$fit + std * newdat$se.fit)
#' newdat$fit <- fit$family$linkinv(newdat$fit)  # Rescale to 0-1
#'
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' lines(pmat ~ L, newdat, col=8, lty=3)
#' polygon(
#'   x = c(newdat$L, rev(newdat$L)),
#'   y = c(newdat$ymax, rev(newdat$ymin)),
#'   col = adjustcolor(2, alpha.f = 0.3),
#'   border = adjustcolor(2, alpha.f = 0.3)
#' )
#' lines(fit ~ L, newdat, col=2)
#'
#' lrPerc <- function(alpha, beta, p) (log(p/(1-p))-alpha)/beta
#' ( L50 <- lrPerc(alpha=coef(fit)[1], beta=coef(fit)[2], p=0.5) )
#' lines(x=c(L50,L50,0), y=c(-100,0.5,0.5), lty=2, col=2)
#' text(x=L50, y=0.5, labels = paste0("L50 = ", round(L50,2)), pos=4, col=2 )
#'
#'
#' }
#'
#'



virtualPop <- function(
tincr = 1/12,
K.mu = 0.5, K.cv = 0.1,
Linf.mu = 80, Linf.cv = 0.1,
ts = 0, C = 0.85,
LWa = 0.01, LWb = 3,
Lmat = 40, wmat = 8,
rmaxBH = 10000,
betaBH = 1, srr.cv = 0.1,
repro_wt = c(0,0,0,1,0,0,0,0,0,0,0,0),
M = 0.7,
Etf = 500,
qtf = 0.001,
harvest_rate = NaN,
gear_types = "trawl",   # alternative: "gillnet"
sel_list = list(mesh_size=100, mesh_size1=60,select_dist="lognormal",select_p1=3, select_p2=0.5),  # parameters adapted from the tilapia data set (increased spread)
L50 = 0.25*Linf.mu,
wqs = L50*0.2,
bin.size = 1,
timemin = 0, timemax = 20, timemin.date = as.Date("1980-01-01"),
N0 = 10000,
fished_t = seq(17,20,tincr),
lfqFrac = 1,
progressBar = TRUE
){

    ## Fishing mortality - effort - catchability

    ## if E == single value, assuming one fleet and same effort for all fished years
    if(length(as.numeric(Etf))==1){
        Emat <- as.matrix(rep(Etf,length(fished_t)))
    }
    ## if E == matrix, rows = years and columns = fleets
    if(class(Etf) == "matrix"){
        Emat <- Etf
    }else if(length(Etf)>1){
        Emat <- as.matrix(Etf)
    }

    ## adapt q to Emat
    if(length(qtf)==1){
        qmat <- matrix(qtf, ncol=dim(Emat)[2], nrow=dim(Emat)[1])
    }
    if(class(qtf) == "matrix"){
        if(dim(qtf)[1] != dim(Emat)[1]){
            qmat <- matrix(rep(qft[1,], dim(Emat)[1]),ncol=dim(qft)[2],byrow=TRUE)
        }else{
            qmat <- qtf
        }
    }else if(length(qtf)>1){
        qmat <- as.matrix(qtf)
    }

    ## If no harvest_rate provided assuming that effort * catchability = fishing mortality
    if(is.na(harvest_rate) | is.nan(harvest_rate)){
        harvest_rate <- Emat * qmat
    }

    selfunc <- function(Lt, fleetNo){
        if(is.na(fleetNo)){
            gear_typesX <- gear_types
            L50X <- L50
            wqsX <- wqs
            sel_listX <- sel_list
        }else{
            gear_typesX <- gear_types[fleetNo]
        }
        switch(gear_typesX,
               trawl ={
                   if(!is.na(fleetNo)){
                       L50X <- L50[fleetNo]
                       wqsX <- wqs[fleetNo]
                   }
                   pSel <- logisticSelect(Lt=Lt, L50=L50X, wqs=wqsX)},
               gillnet={
                   if(!is.na(fleetNo)){
                       sel_listX <- sel_list[[fleet_No]]
                   }
                   pSel <- do.call(fishdynr::gillnet, c(list(Lt=Lt),sel_listX))
               },
               stop(paste("\n",gear_typesX,"not recognized, possible options are: \n","trawl \n","gillnet \n")))
        return(pSel)
    }

    ## ## if multiple fleets target the same stock, the harvest rate of each fleet is scaled according to the combined harvest rate - this only works if all fleets would have the same gear!
    ## if(class(harvest_rate) == "matrix"){
    ##     multimat <- harvest_rate / rowSums(harvest_rate)
    ##     harvest_rate <- rowSums(harvest_rate * multimat)
    ##    }

# times
timeseq = seq(from=timemin, to=timemax, by=tincr)
if(!zapsmall(1/tincr) == length(repro_wt)) stop("length of repro_wt must equal the number of tincr in one year")
repro_wt <- repro_wt/sum(repro_wt)
repro_t <- rep(repro_wt, length=length(timeseq))
# repro_t <- seq(timemin+repro_toy, timemax+repro_toy, by=1)

# make empty lfq object
lfq <- vector(mode="list", length(timeseq))
names(lfq) <- timeseq

indsSamp <- vector(mode="list", length(timeseq))
names(indsSamp) <- timeseq


# Estimate tmaxrecr
tmaxrecr <- (which.max(repro_wt)-1)*tincr

# mean phiprime
phiprime.mu = log10(K.mu) + 2*log10(Linf.mu)



# required functions ------------------------------------------------------
date2yeardec <- function(date){as.POSIXlt(date)$year+1900 + (as.POSIXlt(date)$yday)/365}
yeardec2date <- function(yeardec){as.Date(strptime(paste(yeardec%/%1, ceiling(yeardec%%1*365+1), sep="-"), format="%Y-%j"))}

make.inds <- function(
	id=NaN, A = 0, L = 0, W=NaN, mat=0,
	K = K.mu, Winf=NaN, Linf=NaN, phiprime=NaN,
	F=NaN, Z=NaN,
	Fd=0, alive=1
){
  inds <- data.frame(
    id = id,
    A = A,
    L = L,
    W = W,
    Lmat=NaN,
    mat = mat,
    K = K,
    Linf = Linf,
    Winf = Winf,
    phiprime = phiprime,
    F = F,
    Z = Z,
    Fd = Fd,
    alive = alive
  )
  lastID <<- max(inds$id)
  return(inds)
}

express.inds <- function(inds, seed){
  set.seed(seed)
  inds$Linf <- Linf.mu * rlnorm(nrow(inds), 0, Linf.cv)
  inds$Winf <- LWa*inds$Linf^LWb
  # inds$K <- 10^(phiprime.mu - 2*log10(inds$Linf)) * rlnorm(nrow(inds), 0, K.cv)
  seed1 <- seed + 1
  set.seed(seed1)
  inds$K <- K.mu * rlnorm(nrow(inds), 0, K.cv)
  inds$W <- LWa*inds$L^LWb
  inds$phiprime <- log10(inds$K) + 2*log10(inds$Linf)
  seed2 <- seed + 2
  set.seed(seed2)
  inds$Lmat <- rnorm(nrow(inds), mean=Lmat, sd=wmat/diff(qnorm(c(0.25, 0.75))))
	return(inds)
}

grow.inds <- function(inds){
	# grow
  L2 <- dt_growth_soVB(Linf = inds$Linf, K = inds$K, ts = ts, C = C, L1 = inds$L, t1 = tj-tincr, t2 = tj)
  # update length and weight
	inds$L <- L2
	inds$W <- LWa*inds$L^LWb
	# age inds
	inds$A <- inds$A + tincr
	return(inds)
}

mature.inds <- function(inds){
	# p <- pmat_w(inds$L, Lmat, wmat) # probability of being mature at length L
	# p1t <- 1-((1-p)^tincr)
	# inds$mat <- ifelse(runif(nrow(inds)) < p1t | inds$mat == 1, 1, 0)
  inds$mat <- ifelse((inds$L > inds$Lmat | inds$mat == 1), 1, 0)
	return(inds)
}

death.inds <- function(inds, f0 = FALSE){
        ## multiple fleets
        if(class(harvest_rate)=="matrix"){
            if(dim(harvest_rate)[2]>=2){
                pSel <- matrix(NaN, ncol=dim(harvest_rate)[2],nrow=dim(inds)[1])
                for(seli in 1:(dim(harvest_rate)[2])){
                    pSel[,seli] <- selfunc(Lt = inds$L, fleetNo = seli)
                }
                ## effective fishing mortality (in relation to selectivity) - per fleet with mutliple fleets
                Feff <- pSel * Fmax
                ## single fishing mortality value (per year) scaled according to F of each fleet
                ## this calculation only works if there is fishery (Fmax in denominator not allowed to be 0, otherwise F = NaN and then Z = NaN), thus:
                if(all(Fmax == 0)){
                    inds$F <- 0
                }else{
                    inds$F <- as.numeric(rowSums(Feff * Fmax) / sum(Fmax))
                }
            }else{
            ## single fleet
                pSel <- selfunc(Lt = inds$L, fleetNo = NA)
                inds$F <- as.numeric(pSel * Fmax)
            }
        }else{
            ## single fleet
            pSel <- selfunc(Lt = inds$L, fleetNo = NA)
            inds$F <- as.numeric(pSel * Fmax)
        }
        inds$Z <- M + inds$F
        if(f0) inds$Z <- M
	pDeath <- 1 - exp(-inds$Z*tincr)
	dead <- which(runif(nrow(inds)) < pDeath)
	# determine if natural or fished
	if(length(dead) > 0){
	  inds$alive[dead] <- 0
	  tmp <- cbind(inds$F[dead], inds$Z[dead])
	  # Fd=1 for fished individuals; Fd=0, for those that died naturally
	  Fd <- apply(tmp, 1, FUN=function(x){sample(c(0,1), size=1, prob=c(M/x[2], x[1]/x[2]) )})
  	inds$Fd[dead] <- Fd
    rm(tmp)
	}
	return(inds)
}

remove.inds <- function(inds){
  dead <- which(inds$alive == 0)
  if(length(dead)>0) {inds <- inds[-dead,]}
  return(inds)
}

reproduce.inds <- function(inds, seed){
	## reproduction can only occur of population contains >1 mature individual
    if(repro > 0 & sum(inds$mat) > 0){
        ## calc. SSB
        SSB <- sum(inds$W*inds$mat)
        n.recruits <- ceiling(srrBH(rmaxBH, betaBH, SSB) * repro)
        ## add noise to recruitment process
        seed3 <- seed + 3
        set.seed(seed3)
        n.recruits <- n.recruits * rlnorm(1, 0, sdlog = srr.cv)
        ## make recruits
        offspring <- make.inds(
            id = seq(lastID+1, length.out=n.recruits)
        )
        ## express genes in recruits
        offspring <- express.inds(offspring, seed = seed + 4)
        ##combine all individuals
        inds <- rbind(inds, offspring)
    }
    return(inds)
}

record.inds <- function(inds, ids=1:10, rec=NULL){
	if(is.null(rec)) {
		rec <- vector(mode="list", length(ids))
		names(rec) <- ids
		inds <- inds
	} else {
		ids <- as.numeric(names(rec))
	}
	if(length(rec) > 0) {
		inds.rows.rec <- which(!is.na(match(inds$id, ids)))
		if(length(inds.rows.rec) > 0){
			for(ii in inds.rows.rec){
				match.id <- match(inds$id[ii], ids)
				if(is.null(rec[[match.id]])) {
					rec[[match.id]] <- inds[ii,]
				} else {
					rec[[match.id]] <- rbind(rec[[match.id]], inds[ii,])
				}
			}
		}
	}
	return(rec)
}

spm <- function(x, B0){
  K = x[1]
  r = x[2]
  n = x[3]
  Bthat <- rep(NA, length(resf0$pop$B))
  Bthat[1] <- B0
  for(i in 2:length(Bthat)){
    Bthat[i] <- Bthat[i-1] + (r / (n - 1)) * Bthat[i-1] * (1 - (Bthat[i-1] / K)^(n-1))
  }
  sum((resf0$pop$B - Bthat)^2)
}

spm2 <- function(pars){
  B0 <- pars[1]
  K <- pars[2]
  r <- pars[3]
  n <- pars[4]
  Bthat <- rep(NA, length(resf0$pop$B))
  Bthat[1] <- B0
  for(i in 2:length(Bthat)){
    Bthat[i] <- Bthat[i-1] + (r / (n - 1)) * Bthat[i-1] * (1 - (Bthat[i-1] / K)^(n-1))
  }
  Bthat
}


# run model ---------------------------------------------------------------

# Initial population
lastID <- 0
inds <- make.inds(
  id=seq(N0)
)
inds <- express.inds(inds = inds, seed = 123)

## results object
res <- list()
res$pop <- list(
  dates = yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin) ),
	N = NaN*timeseq,
	B = NaN*timeseq,
	SSB = NaN*timeseq
)

## For simulation of unfished population
## Initial population
lastID <- 0
indsf0 <- make.inds(
  id=seq(N0)
)
indsf0 <- express.inds(inds = indsf0, seed = 123)
## results object 
resf0 <- list()
resf0$pop <- list(
  dates = yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin) ),
  N = NaN*timeseq,
  B = NaN*timeseq,
  SSB = NaN*timeseq
)

## simulation
if(progressBar) pb <- txtProgressBar(min=1, max=length(timeseq), style=3)
for(j in seq(timeseq)){
  tj <- timeseq[j]
  
  ## harvest rate applied? lfq sampled?
  if(is.na(fished_t[1]) | is.nan(fished_t[1])){ ## before: length(fished_t) == 0 : as I see it fished_t never has length 0, even if set ot NA or NaN, it woudl have length 1
    Fmax <- 0
    lfqSamp <- 0
  } else if(min(sqrt((tj-fished_t)^2)) < 1e-8){
       ## time index for fished_t
          tfish <- which.min(abs(fished_t - tj))
          ## provide yearly Fmax value (per fleet if multiple fleets simulated)
          if(class(harvest_rate) == "matrix"){
              Fmax <- harvest_rate[tfish,]
          }else if(length(harvest_rate)>1){
             Fmax <- harvest_rate[tfish]
          }else{
              Fmax <- harvest_rate
          }
          lfqSamp <- 1
          } else {
      Fmax <- 0
      lfqSamp <- 0
    }

  repro <- repro_t[j]

	# population processes
	inds <- grow.inds(inds)
	inds <- mature.inds(inds)
	inds <- reproduce.inds(inds = inds, seed = j)
	inds <- death.inds(inds)
	## sample lfq data

	if(lfqSamp){
	  samp <- try( sample(seq(inds$L), ceiling(sum(inds$Fd)*lfqFrac), prob = inds$Fd), silent = TRUE)

    ## tmp <- try( sample(inds$L, ceiling(sum(inds$Fd)*lfqFrac), prob = inds$Fd), silent = TRUE)
    if(class(samp) != "try-error"){
      lfq[[j]] <- inds$L[samp]
      indsSamp[[j]] <- inds[samp,]
    }
    rm(samp)
  }
	inds <- remove.inds(inds)

	# update results
	res$pop$N[j] <- nrow(inds)
	res$pop$B[j] <- sum(inds$W)
	res$pop$SSB[j] <- sum(inds$W*inds$mat)
	
	## simulate unfished population for K, r and SSB_F=0
	# population processes
	indsf0 <- grow.inds(indsf0)
	indsf0 <- mature.inds(indsf0)
	indsf0 <- reproduce.inds(inds = indsf0, seed = j)
	indsf0 <- death.inds(indsf0, f0 = TRUE)
	indsf0 <- remove.inds(indsf0)
	
	# update results
	resf0$pop$N[j] <- nrow(indsf0)
	resf0$pop$B[j] <- sum(indsf0$W)
	resf0$pop$SSB[j] <- sum(indsf0$W * indsf0$mat)

	## update progressbar
	if(progressBar) setTxtProgressBar(pb, j)
}
if(progressBar) close(pb)


## Estimate carrying capacity
## only if years without fishing are simulated (subsequently to the first 10 years)
    # if(length(fished_t) < length(timeseq)){
    #     startyear  <- as.POSIXlt(timemin.date)
    #     startyear$year <- startyear$year + 10
    #     year10 <- as.Date(startyear)
    #     cutoff <- which.min(abs(yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin)) - year10))
    #     cc_years <- seq(timeseq)[-c((1:cutoff),which(round(timeseq,5) %in% round(fished_t, 5)))]
    #     if(length(cc_years) > 3){
    #         mod <- lm(res$pop$B[cc_years] ~ 1)
    #         res$pop$K <- as.numeric(coefficients(mod))
    #     }
    # }
## Alternative way build into loop above
startyear  <- as.POSIXlt(timemin.date)
startyear$year <- startyear$year + 10
year10 <- as.Date(startyear)
cutoff <- which.min(abs(yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin)) - year10))
cc_years <- seq(timeseq)[-(1:cutoff)]
if(length(cc_years) > 3){
      mod <- lm(resf0$pop$B[cc_years] ~ 1)
      resf0$pop$K <- as.numeric(coefficients(mod))
}
if(length(cc_years) > 3){
  mod <- lm(resf0$pop$SSB[cc_years] ~ 1)
  resf0$pop$SSBf0 <- as.numeric(coefficients(mod))
}

## estimate K, r, n 
resSPM <- optim(par = c(resf0$pop$K, 0.5, 2), fn = spm, B0 = resf0$pop$B[1])
Kest <- resSPM$par[1]
rest <- resSPM$par[2]
nest <- resSPM$par[3]
gammal <- nest^(nest/(nest-1))/(nest-1)
m <- rest * Kest / (nest^(nest/(nest-1)))
## Deterministic reference levels
if(nest == 1){
  ## Fox reference levels
  Bdmsy <- Kest/exp(1)
  msyd <- rest * Kest / exp(1)
  Fdmsy <- rest
}else{
  ## Pella and Tomlinson reference levels
  Bdmsy <- nest^(1/(1-nest)) * Kest
  msyd <- m
  Fdmsy <- m/Bdmsy
}

## save parameters
resf0$pop$K2 <- Kest
resf0$pop$r <- rest
resf0$pop$n <- nest
resf0$pop$m <- m
resf0$pop$gamma <- gammal
resf0$pop$Bdmsy <- Bdmsy
resf0$pop$msyd <- msyd
resf0$pop$Fdmsy <- Fdmsy

## Production curve
spmPlot <- spm2(c(res$pop$B[1], Kest,  rest, nest))
Prod <- (rest / (nest - 1)) * spmPlot * (1 - (spmPlot / Kest)^(nest-1))

# Export data -------------------------------------------------------------

    ## for simulation of population without exploitation, necessary to make the lfq export optional:
    if(any(!is.na(fished_t[1]) & !is.nan(fished_t[1])) & (lfqFrac != 0 & !is.na(lfqFrac) & !is.nan(lfqFrac))){
        ## Trim and Export 'lfq'
        lfq2 <- lfq[which(sapply(lfq, length) > 0)]
        ## binned version of lfq
        dates <- yeardec2date( date2yeardec(timemin.date) + (as.numeric(names(lfq2)) - timemin) )
        Lran <- range(unlist(lfq2))
        Lran[1] <- floor(Lran[1])
        Lran[2] <- (ceiling(Lran[2])%/%bin.size + ceiling(Lran[2])%%bin.size + 1) * bin.size
        bin.breaks <- seq(Lran[1], Lran[2], by=bin.size)
        bin.mids <- bin.breaks[-length(bin.breaks)] + bin.size/2
        res$lfqbin <- list(
            sample.no = seq(bin.mids),
            midLengths = bin.mids,
            dates = dates,
            catch = sapply(lfq2, FUN = function(x){
                hist(x, breaks=bin.breaks, plot = FALSE, include.lowest = TRUE)$counts
            })
        )
    }


# individuals
indsSamp <- indsSamp[which(sapply(indsSamp, length) > 0)]
res$inds <- indsSamp


# record mean parameters
res$growthpars <- list(
  K = K.mu,
  Linf = Linf.mu,
  C = C,
  ts = ts,
  t_anchor = weighted.mean(date2yeardec(as.Date(paste("2015",which(repro_wt != 0),"15",sep="-"))) %% 1,
              w =repro_wt[which(repro_wt != 0)]),  ## weighted mean of t_anchor
  phiprime = phiprime.mu,
  tmaxrecr = tmaxrecr
)

    ## fisheries dependent information
    ## if fisheries are simulated
    if(any(!is.na(fished_t) & !is.nan(fished_t))){
        res$fisheries <- list(
            fished_years = yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin) )[fished_t],
            E = Emat,
            q = qmat,
            F = harvest_rate
        )
    }

return(res)

if(plot){
  opar <- par(mfrow=c(3,2))
  ## Numbers
  with(res$pop, plot(dates, N, type='l', lwd=2, 
                     xlab="",ylab="Numbers",
                     main = "Population trajectory",
                     ylim=c(0,max(resf0$pop$N,na.rm=TRUE))))
  with(resf0$pop, lines(dates, N, col=4, lwd=2))
  points(res$pop$dates[1], N0,pch=4, lwd=2, col='darkred')
  ## Biomass  + SSB
  with(res$pop, plot(dates, B, type='l', lwd=2, 
                     xlab="",ylab="Biomass",
                     main = "Biomass trajectory",
                     ylim=c(0,max(resf0$pop$B,na.rm=TRUE))))
  with(res$pop, lines(dates, SSB, lwd=2, lty=3))
  with(resf0$pop, lines(dates, B, col=4, lwd=2))
  with(resf0$pop, lines(dates, SSB, col=4, lwd=2, lty=3))
  abline(h = resf0$pop$K, col='darkred', lwd=2)
  abline(h = resf0$pop$SSBf0, col='darkred', lwd=2,lty=3)
  ## SPM
  with(resf0$pop, plot(dates, B, type='b', lwd=2,
                     xlab="",ylab="Biomass",col=4,
                     pch=16,
                     main = "Surplus production model",
                     ylim=c(0,max(resf0$pop$B,na.rm=TRUE))))
  abline(h = resf0$pop$K2,lwd=2, lty=3, col = 'darkred')
  with(resf0$pop, lines(dates, spmPlot, lwd=2, lty=1, col = 'darkred'))
  ## Production curve
  plot(spmPlot, Prod, type='l', lwd=2, col = 4, 
       main = "Production curve",
       xlab="Biomass",ylab="Surplus production", ylim = c(0,max(Prod,na.rm=TRUE)*1.1))
  segments(x0 = 0, y0 = msyd, x1 = Bdmsy, y1 = msyd, col = "darkred", lty=3, lwd=2)
  segments(x0 = Bdmsy, y0 = 0, x1 = Bdmsy, y1 = msyd, col = "darkred", lty=3, lwd=2)
  
  
  
  
  ## Legend
  plot.new()
  legend("center", legend=c("fishing", "no fishing", "population parameters"),
         col=c("black","blue", "darkred"), lwd=2, lty=1, bty='n', 
         y.intersp = 0.4, x.intersp = 0.3, seg.len = 0.2)
  
  par(opar)
}
} # end of function


res$pop$B == res$pop$SSB
