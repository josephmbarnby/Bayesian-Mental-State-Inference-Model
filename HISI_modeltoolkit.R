#######  MODELS & MAP WRAPPERS ###############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(matlab)
library(R.matlab)
library(tidyverse)
library(tidyquant)
library(ggcorrplot)

### Dictator_Reversal_Functions.R ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

infHISIll_20d <- function(parms,d,details=0,plots=0,sim_only=1,tn = 10,phase = 2) {

  d <- as.matrix(d)

  if(sim_only == 1){
  tn    = tn;  # Trial number per phase
  } else if (sim_only == 0){
  tn    = length(d[,1])/phase
  }

  phase = phase;  # number or partners or phases

  Nb = 9;     # number of bins for each personality dim.
  Nf = 9;     # 'full' number of action bins (within which density is flat) for detailed sim work
  Na = 2;     # actual number of actions from a partner

  # i.e. proportions that may be received by the Other.
  # In general will be set to Nb.

  #   Prior beliefs of the pt ** about the population **
  #   Prior beliefs about greed on its own, and malice on its own:
  PSI0 = noisyBino(parms[3], parms[4],Nb); PHI0 = noisyBino(parms[1],parms[2],Nb);
  # In this version, these baseline beliefs are considered as independent,
  # so a simple matrix multiplication gives the joint:
  PSIHI0 = PSI0 %*% t(PHI0);

  if (plots){
    opar=par()
    anames = c('selfish','fair')
  }

  # Now to formulate the policies. There is one individual parameter
  # here, which determines the uncertainty with which say a high-HI person will
  # indeed choose the action most fitting to their type (i.e., keep everything),
  # or will show choices more variable over options:
  upi = parms[5];     # convenience copy of policy variability (inverse precision)
  w0  = parms[6];     # This and two next will help map HI, SI onto policy (zero or half return)
  whi = parms[7];
  wsi = parms[8];

  if ((whi < 0) || (wsi < 0)){
    print('Params c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi,etaHI (or single one), (optionally etaSI)) :')
    print(parms)
    error("whi,wsi must all be non-negative")
  }

  if (length(parms) == 8){
  etaHI = 1
  etaSI = 1
  } else if (length(parms) == 9){
  eta = parms[9];
  } else if (length(parms) == 10){
  etaHI = parms[9];
  etaSI = parms[10]
  }

  err =0.02/(Nb*Nb) # arbitrary lapse-rate-like
  #}
  # param; Note scaling by the number of attribution states considered.

  # Set up the map between 'attributes' and actions :
  pi  = array(NA,c(Nb,Nb,Na));

  offs = (Nb+1)/2
  for (SI in 1:Nb){
    for (HI in 1:Nb){
      pi[SI,HI,1]  = invlogit(w0  + (wsi*(SI-offs)) + (whi*(HI-offs)))  # prob. of ufair offer goes up w. HI, SI
      # the offs is to center at (usually) 5
      pi[SI,HI,2]  = 1 - pi[SI,HI,1]
    }
  }

  ### Run the inference

  # Here the posterior of one trial will form the prior for the next trial,
  # staring from the population prior beliefs PSIHI0.
  #

  sll = 0;
  pri0 = PSIHI0; # prior at the very start of encounter with 'partner'.
  # Construct an output/results object, which is just the sum log lik
  # if details==0. If not, record trial-by-trial data, attributions,
  # prob. of data given params (likelihood), log-lik, most likely attrib given
  # the parameters and data, and also a set of simulated data (incl. decision variability)

  if (details ){
    llout = list();

    llout[[1]]=0;

    # [[2]] is the parameter vector
    hd <- c('pHI0','uHI0','pSI0','uSI0','upi','w0', 'wHI', 'wSI', 'etaHI','etaSI', 'err')
    llout[[2]] = c(parms[1:10],err);  names(llout[[2]]) <- hd;

    # [[3]] is the detailed evolution, evo:
    llout[[3]] = matrix(NA,tn*phase+1,10);
    llout[[3]][,1] = c(0,1:(tn*phase))
    colnames(llout[[3]]) = c('trial','ret','HI','SI','lik','ll','HImode','SImode','HIsim','SIsim')
    llout[[3]][2:(1+tn*phase),2:4] <- as.matrix(d)

    # [[4]] will be a big 9 x 9 x 21 array with detailed policies
    # Hypothetical (attribn. reporting) policy before any data seen:
    pol = pri0^(1/upi); pol = pol / sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    llout[[4]] <- array(dim=c(dim(pri0),(1+tn*phase)))
    llout[[4]][,,1] <- pri0;
    llout[[5]] <- list();
    llout[[6]] <- list();
    names(llout) <- c('sll','par', 'evo','policy','attr2pol', 'plots')
    }

    post <- pri0; # this is the belief that will be updated with each trial

    # rows to be processed and convenience copies of other's actions,
    # malice and greed attributions:
    ro = 1:(phase*tn); # rows of data matrix
    as = d[ro,1];  aind = round((Na-1)*as+1)
    hi = d[ro,2];  hind = round((Nb-1)*hi+1)
    si = d[ro,3];  sind = round((Nb-1)*si+1)

    for (t in 1:(tn*phase)){  # loop

      #if(plots){
      #  heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
      #          margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
      #          main = paste('\n lnPost. at start of trial ',t),xlab='HI',ylab='SI')
      #}


      if (t == (tn+1) & length(parms)==9){

        post = (pri0 * (1-eta)) + (post * eta);

      } else if (t == (tn+1) & length(parms)==10) {

        PSIPost <- rowSums(post);
        PHIPost <- colSums(post);
        post    <- ((PSI0 * (1-etaSI)) + (etaSI * PSIPost)) %*% t((PHI0 * (1-etaHI)) + (etaHI * PHIPost))        ;

      }

      pri = post;              # new prior is last posterior

      # In the next line, the pt. uses the pi entry as likelihood, pri as prior,
      # over the character of the partner. This post is their post. beliefs
      post = pi[,,aind[t]] * pri
      post = post / sum(as.vector(post))  # Bayes

      # Now the probability of the response, incl. the lapse-rate-like err:
      pol = post^(1/upi);
      pol = pol/sum(as.vector(pol)); #renormalise

      pol = (pol+err)/(1+err*length(pol))
      lik = pol[sind[t],hind[t]];
      sll = sll + log(lik);         # accumulate sum log lik

      if (details ){

        llout$evo[(t+1),'lik'] <- lik
        llout$evo[(t+1),'ll'] <- log(lik)
        # find mode of pol
        c = max.col(pol);  # this finds the max. col. pos. for each row
        m=c(); for (r in 1:Nb){m[r]=pol[r,c[r]]};  # retrieve the values ...
        r=which.max(m);  # ... so was to find the row with the mode of the distro.
        llout$evo[(t+1),c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
        # Now sample randomly
        hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
        sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the phase from the corresponding conditional.
        # debug: barplot(colSums(pol));       barplot(pol[,hisim]);
        llout$evo[(t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
        llout$policy[,,(t+1)] <- pol

      }
    }

    if(plots){
      library(tidyquant)
      library(ggplot2)
      library(patchwork)
      initial_policy <- pi[,,1] %>%
        as.data.frame() %>%
        mutate(SI = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85)) %>%
        pivot_longer(1:9, 'HI', values_to = 'Probability') %>%
        ggplot(aes(HI, SI, fill = Probability)) +
        geom_tile() +
        scale_fill_gradient2(low = '#98C1D9', mid = 'white', high = '#C81D25', midpoint = 0.5,
                             labels = seq(0, 1, 0.2), breaks = seq(0, 1, 0.2))+
        scale_x_discrete(labels = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85))+
        theme_tq()
      plotobject <- llout$evo %>%
        as.data.frame() %>%
        na.omit()
      if(sim_only == 1){
      plotobject2 <- plotobject  %>%
        pivot_longer(HIsim:SIsim, 'Attribution', values_to = 'values') %>%
        mutate(ret = ifelse(ret == 0.5, 1, 0))
      trial_wise <- ggplot( plotobject2, aes(trial, values, color = Attribution))+
          geom_line()+
          geom_point(aes(trial, ret), color = ifelse( plotobject2$ret == 1, 'black', 'dark red'))+
          scale_color_brewer(palette = 'Set1')+
          labs(x = 'Trial', y = 'Simulated Observation')+
          theme_tq()
      }else if (sim_only == 0){
      trial_wise <- plotobject  %>%
        pivot_longer(HIsim:SIsim, 'Attribution', values_to = 'values') %>%
        pivot_longer(HI:SI, 'AttributionREAL', values_to = 'valuesREAL') %>%
        mutate(ret = ifelse(ret == 0.5, 1, 0)) %>%
        ggplot(aes(trial, values, color = Attribution))+
          geom_line(linetype = 1)+
          geom_line(aes(trial, valuesREAL, color = AttributionREAL), linetype = 2)+
          geom_point(aes(trial, ret), color = ifelse(ret == 1, 'black', 'dark red'))+
          scale_color_brewer(palette = 'Set1')+
          labs(x = 'Trial', y = 'Simulated Observation')+
          theme_tq()
      }
      fitness <- plotobject %>%
        ggplot(aes(ll))+
        geom_density()+
        geom_histogram(alpha = 0.5, binwidth = 0.2)+
        geom_vline(xintercept = -4.17)+
        labs(x = 'Trial-wise LL', y = 'Density')+
        theme_tq()
      plotter <- (trial_wise | (initial_policy/fitness)) + plot_annotation(tag_levels = 'A')
      llout[[6]] <- plotter
    }

    #if (plots) {
    #  heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
    #          margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
    #          main = paste('\n lnPost. after triall ',t),xlab='HI',ylab='SI')
    #  par(opar)  # restore state of graphics where we found it.
    #}

  if (details ){llout$sll <- sll} else {llout <- sll}
  return(llout)

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.a. Wrapper for function wth attribution --> policy params

# params c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi,etaHI (or single one), (optionally etaSI) ...
msLPhisi_20d <- function(ParM, datAr, scbeta0=-1,details=0){

  parM <- as.vector(ParM) # in case it's inputted in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.01,1.01,0, 1,     # 1  for pHI0
                        1, 3.6, 0, 25,    # 2  for uHI0
                        1.01,1.01,0, 1,     # 3  for pSI0
                        1, 3.6, 0, 25,    # 4  for uSI0
                        1.01,1.01, 0,1,    # 5  for upi
                        3.6, 3.6, -25, 25,  # 6  for w0
                        1.01,1.01, 0, 1,    # 7  for whi
                        1.01,1.01, 0, 1,    # 8  for wlo
                        1.01,1.01,0, 1,     # 9  for eta/etaHI
                        1.01,1.01,0, 1      # 10 for etaSI
                        ),
                        nrow=4, ncol=10)
    if(parn==10){
      scbeta0 = scbeta0
    }else if (parn == 9) {
      scbeta0 = scbeta0[,1:9]
    } else {
      scbeta0 = scbeta0[,1:8]
    }

    if(details){
      colnames(scbeta0) <-   c('pHI0','uHI0','pSI0','uSI0','upi','w0',   'whi',  'wsi','etaHI', 'etaSI')
      rownames(scbeta0) <-   c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have say 24 elements or more!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - infHISIll_20d(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
      res[[1]] <- Inf; res[[5]] <- NA;
    } else {
      res[[5]] <- infHISIll_20d(ParM,datAr);
      res[[1]] <- mSLPrior - res[[5]];
    }
    names(res) <- c('sLP','scbeta0','par','dat','sLL')
    return(res)
  }


} # end of msLPhisi_20d

# Generate simulations ----------------------------------------------------

simulatedata_HISI <- function(x = 'none',
                              values = 10,
                              samples = 100,
                              plot = 1,
                              trials = 10,
                              partners = 2,
                              partner_type = 'fair',
                              prob = 0.2,
                              pHI0 = 0.5, uHI0 = 2, pSI0 = 0.5, uSI0 = 2, w0 = 0, wHI = 0.1, wSI = 0.1, upi = 2, eta = 0.5) {

  if(!x %in% c('pHI0', 'uHI0', 'pSI0', 'uSI0', 'upi', 'upi01', 'w0', 'wHI', 'wSI', 'eta', 'none', 'wHI01', 'wSI01', 'prob')){
    stop("Incorrect input! Must be one of: 'pHI0', 'uHI0', 'pSI0', 'uSI0', 'upi', 'upi01', 'w0', 'wHI', 'wSI', 'eta', 'none', 'wHI01', 'wSI01', 'prob'")
  }

  if(pHI0 > 1) stop('pHI0 must be < 1')
  if(pSI0 > 1) stop('pSI0 must be < 1')
  if(eta  > 1) stop('eta must be < 1')


  if(x == 'none'){values = 1}

  probA = prob
  probB = 1-probA
  fair_d   = sample(c(0.5, 0), trials, prob = c(probB, probA), replace = T)
  unfair_d = sample(c(0.5, 0), trials, prob = c(probA, probB), replace = T)
  partial_d= sample(c(0.5, 0), trials, prob = c(0.5, 0.5), replace = T)

  d <- matrix(0, nrow = trials * partners, ncol = 3)

  test_upi    <- list()
  loop_k      <- list()
  loop_order  <- list()
  list_produce<- list()

  for (order in 1:partners){
    for (k in 1:values){
      for (i in 1:samples){

      if        (order == 1 & partners == 2) {d[,1] = c(unfair_d, fair_d)
      #} else if (order == 1 & partners == 2 & partner_type == 'fair') {d[,1] = c(fair_d, unfair_d)
      } else if (order == 2 & partners == 2) {d[,1] = c(fair_d, unfair_d)
      } else if (order == 1 & partners == 3) {d[,1] = c(unfair_d, fair_d, partial_d)
      } else if (order == 2 & partners == 3) {d[,1] = c(fair_d, partial_d, unfair_d)
      } else if (order == 3 & partners == 3) {d[,1] = c(partial_d, unfair_d, fair_d)
      } else if (order == 1 & partners == 1 & partner_type == 'fair'){d[,1] = c(fair_d)
      } else if (order == 1 & partners == 1 & partner_type == 'unfair'){d[,1] = c(unfair_d)
      } else if (order == 1 & partners == 1 & partner_type == 'random'){d[,1] = c(partial_d)
      }

      if(x == 'pHI0') {meanobject = infHISIll_20d(c(k/values, uHI0, pSI0    , uSI0, upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'uHI0') {meanobject = infHISIll_20d(c(pHI0    , k   , pSI0    , uSI0, upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'pSI0') {meanobject = infHISIll_20d(c(pHI0    , uHI0, k/values, uSI0, upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'uSI0') {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , k   , upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'upi')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, k       , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'upi01'){meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, k/values, w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'w0')   {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , k , wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'wHI')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, k       , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'wSI')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, wHI     , k       , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'wHI01'){meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, k/values, wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'wSI01'){meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, wHI     , k/values, eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'eta')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, wHI     , wSI     , k/values), d, details = 1, tn = trials, phase = partners)}
      if(x == 'none') {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, wHI     , wSI     , eta),      d, details = 1, tn = trials, phase = partners)}
      if(x == 'prob' & partners == 1) {
        probA = k/10
        probB = 1-probA
        d[,1] <- sample(c(0.5, 0), trials, prob = c(probB, probA), replace = T)
        meanobject = infHISIll_20d(c(pHI0, uHI0, pSI0, uSI0, upi, w0, wHI, wSI, eta), d, details = 1, tn = trials, phase = partners)
        }

      loop_k[[i]]     <- meanobject$evo %>%
                            as.data.frame() %>%
                            na.omit() %>%
                            mutate(loop = ifelse(x %in% c('pHI0', 'pSI0', 'wHI01', 'wSI01', 'upi01', 'eta', 'prob'), k/values, k),
                                   id = i,
                                   Order = order)

      }

    loop_order[[k]] <- do.call(rbind,loop_k)

    if(samples>1){
    cat(paste('\n Simulated ', k*samples, 'iterations out of ',values*samples,'| Order = ', order, '| Parameter = ', x))
    }
    }

  test_upi[[order]] <- do.call(rbind, loop_order)

  }#end

  list_produce[[1]] <- do.call(rbind, test_upi)

    if(plot){
    list_produce[[2]] <- list_produce[[1]]  %>%
      pivot_longer(HIsim:SIsim, 'Attribution', values_to = 'values') %>%
      mutate(Attribution = ifelse(Attribution == 'HIsim', 'HI Simulated', 'SI Simulated'),
             Order = ifelse(Order == 1 & partners == 2, 'Unfair | Fair',
                            ifelse(Order == 2 & partners == 2, 'Fair | Unfair',
                                   ifelse(Order == 1 & partners == 3, 'Unfair | Fair | Partial',
                                          ifelse(Order == 2 & partners == 3, 'Fair | Partial | Unfair',
                                                 ifelse(Order == 3 & partners == 3, 'Partial | Unfair | Fair', NA))))),
             ret = ifelse(ret == 0.5, 1, 0))

     list_produce[[2]] <- ggplot(list_produce[[2]], aes(trial, values, color = loop, group = loop, linetype= Attribution))+
        stat_summary(geom = 'line')+
        scale_color_gradient(low = '#BFD7EA', high = '#0B3954', name = x)+
        facet_wrap(Attribution~Order, scales = 'free_y', nrow = 2, ncol = partners)+
        coord_cartesian(ylim = c(0, 1))+
        labs(x = 'Trial', y = 'Simulated Observation')+
        theme_tq()+
        theme(text = element_text(size = 14),
              strip.background.x = element_rect(fill = '#009DDC'),
              legend.key.width = unit(1, 'cm'))

     if(x == 'none'){
     list_produce[[2]] <- list_produce[[2]] +
       scale_color_gradient(guide = 'none') +
       labs(title = expr(paste('pHI'[0], '=', !!pHI0, ' | ',
                      'uHI'[0], '=', !!uHI0,' | ',
                      'pSI'[0], '=', !!pSI0,' | ',
                      'uSI'[0], '=', !!uSI0,' | ',
                      'u'[pi] , '=', !!upi,' | ',
                      'w'[0], '=', !!w0,' | ',
                      'w'[HI],'=', !!wHI,' | ',
                      'w'[SI],'=', !!wSI,' | ',
                      eta,      '=', !!eta))
                      )

     }

     if(partners == 1){
     list_produce[[2]] <- list_produce[[2]] +
       labs(subtitle = paste('Partner = ', partner_type)
                      )+
       facet_wrap(~Attribution, scales = 'free_y', nrow = 2)

     }

     if(partners == 2){

      list_produce[[2]] <- list_produce[[2]]+
        geom_vline(xintercept = c(max(trials)), alpha = 0.5)
     }

     if(partners == 3){

      list_produce[[2]] <- list_produce[[2]]+
        geom_vline(xintercept = c(max(trials), max(trials)*2), alpha = 0.5)
     }

    }

  return(list_produce)

}


# Recovery ----------------------------------------------------------------

RecoverModel_HISI <- function(variations = 100,
                              trials = 10,
                              partners = 2,
                              partner_type = 'fair',
                              corrplot = 1,
                              n_cores = 2){

    registerDoParallel(cores = n_cores)

    recovered <- list()
    permutation <- list()

    for (i in 1:variations){

        pHI0 <- round(mysamp(1, 0.5, 0.2, 0, 1, 1000), 2)
        uHI0 <- round(mysamp(1, 2, 2, 0, 8, 1000), 2)
        pSI0 <- round(mysamp(1, 0.5, 0.2, 0, 1, 1000), 2)
        uSI0 <- round(mysamp(1, 2, 2, 0, 8, 1000), 2)
        upi  <- round(mysamp(1, 2, 2, 0, 1, 1000), 2)
        w0   <- round(mysamp(1, 0, 2, -10, 10, 1000), 2)
        wHI  <- round(mysamp(1, 1, 2, 0, 1, 1000), 2)
        wSI  <- round(mysamp(1, 1, 2, 0, 1, 1000), 2)
        eta  <- round(mysamp(1, 0.5, 0.2, 0, 1, 1000), 2)

        simulated <- simulatedata_HISI(x = 'none', values=1, samples=1, trials=trials, partners=partners, partner_type,
                                       pHI0 = pHI0,
                                       uHI0 = uHI0,
                                       pSI0 = pSI0,
                                       uSI0 = uSI0,
                                       upi = upi,
                                       w0 = w0,
                                       wHI = wHI,
                                       wSI = wSI,
                                       eta = eta,
                                       plot = 0)

        data_sim <- split(simulated[[1]], f = list(simulated[[1]]$id, simulated[[1]]$Order))

        tryP <- c(pHI0, uHI0, pSI0, uSI0, upi, w0, wHI, wSI, eta)

          x <- foreach(pt=1:partners, .combine = rbind) %dopar% {

          fitAttempt <- NA; # clear the decks and attempt the new optim fit
          prm        <- tryP * NA;
          L          <- NA;

          dat          <- data_sim[[pt]][,c(2, 9:10)]
          #opt <- list(NP=500,itermax=100, storepopfrom=1, storepopfreq=1)
          lwr = 0.01
          fitAttempt   <- try(DEoptim(
                                 fn = msLPhisi_20d,
                                 datAr = dat,
                                 scbeta0 = -1,
                                 lower = c(rep(lwr, 5), -10, rep(lwr, 3)),
                                 upper = c(1, 8, 1, 8, 1, 10, 1, 1, 1)))
          if(inherits(fitAttempt, "try-error"))
          {
            #error handling code
          prm          <- tryP # if fit failed, revert to best up to now:
          conv         <- 0
          } else {
          prm          <- as.numeric(fitAttempt$optim$bestmem)
          conv         <- 1
          }

          L <- -msLPhisi_20d(prm,data_sim[[pt]][,c(2, 9:10)],NA)  # plain log-likelihood

          data.frame( ID        =  data_sim[[pt]][1,12],
                      Order     =  data_sim[[pt]][1,13],

                      pHI0_true      =  tryP[1],
                      uHI0_true      =  tryP[2],
                      pSI0_true      =  tryP[3],
                      uSI0_true      =  tryP[4],
                      upi_true       =  tryP[5],
                      w0_true        =  tryP[6],
                      wHI_true       =  tryP[7],
                      wSI_true       =  tryP[8],
                      eta_true       =  tryP[9],

                      pHI0_rec       =  prm[1],
                      uHI0_rec       =  prm[2],
                      pSI0_rec       =  prm[3],
                      uSI0_rec       =  prm[4],
                      upi_rec        =  prm[5],
                      w0_rec         =  prm[6],
                      wHI_rec        =  prm[7],
                      wSI_rec        =  prm[8],
                      eta_rec        =  prm[9],

                      ll         =  L,
                      Variation  =  pt,
                      Convergence=  conv)
          }

    permutation[[i]] <- x
    cat('\n Recovered ', i, ' variations out of ', variations, ' | parms = ', tryP, '||', permutation[[i]]$Convergence)

    }

  recovered[[1]] <- do.call(rbind, permutation)

  if(corrplot == 1){

    recovery_plot <- round(cor(recovered[[1]][,3:20]), 2)
    p.mat         <- ggcorrplot::cor_pmat(recovered[[1]][,3:20])
    recovery_mat_plot  <- ggcorrplot::ggcorrplot(recovery_plot[1:9, 10:18], lab = T, p.mat = p.mat[1:9, 10:18])+
  scale_x_discrete(labels = c(
        pHI0_true      = expression(paste(pHI[0], ' true')),
        uHI0_true      = expression(paste(uHI[0], ' true')),
        pSI0_true      = expression(paste(pSI[0], ' true')),
        uSI0_true      = expression(paste(uSI[0], ' true')),
        upi_true       = expression(paste(u[pi], ' true')),
        w0_true        = expression(paste(w[0], ' true')),
        wHI_true       = expression(paste(w[HI], ' true')),
        wSI_true       = expression(paste(w[SI], ' true')),
        eta_true       = expression(paste(eta, ' true'))
  )) +
  scale_y_discrete(labels = c(
        pHI0_rec      = expression(paste(pHI[0], ' rec')),
        uHI0_rec      = expression(paste(uHI[0], ' rec')),
        pSI0_rec      = expression(paste(pSI[0], ' rec')),
        uSI0_rec      = expression(paste(uSI[0], ' rec')),
        upi_rec       = expression(paste(u[pi], ' rec')),
        w0_rec        = expression(paste(w[0], ' rec')),
        wHI_rec       = expression(paste(w[HI], ' rec')),
        wSI_rec       = expression(paste(w[SI], ' rec')),
        eta_rec       = expression(paste(eta, ' rec'))
  ))

    recovered[[2]] <- recovery_mat_plot

  }

  return(recovered)

}#end


# Fitting -----------------------------------------------------------------

FitModel_HISI <- function(    data = ExampleData,
                              plot = 1,
                              phase = 2,
                              tn = 10,
                              cores = 4){

    registerDoParallel(cores = cores)

    fitted    <- list()
    simulated <- list()

          x <- foreach(pt=1:length(data), .combine = rbind) %dopar% {

          fitAttempt <- NA; # clear the decks and attempt the new optim fit
          prm        <- rep(NA, 9);
          L          <- NA;
          lwr        <- 0.01

          dat         <- data[[pt]][,1:3]
          if(sum(dat[1:tn,1]) > sum(dat[(tn+1):(tn*phase),1])){
            Order_calc <- 'Fair'
          }else {
            Order_calc <- 'Unfair'
          }

          fitAttempt   <- try(DEoptim(
                                     fn = msLPhisi_20d,
                                  datAr = dat,
                                scbeta0 = -1,
                                lower = c(rep(lwr, 5), -10, rep(lwr, 3)),
                                 upper = c(1, 8, 1, 8, 1, 10, 1, 1, 1)))

          if(inherits(fitAttempt, "try-error"))
          {
            #error handling code
          conv = 0
          prm  = rep(NA, 9)# if fit failed, revert to best up to now:
          } else{
          conv = 1
          prm        <- fitAttempt$optim$bestmem;
          }

          L <- -msLPhisi_20d(prm,data[[pt]][,1:3],NA);  # plain log-likelihood

          data.frame( ID        =  ExampleData[[pt]]$ID[1],
                      Order     =  Order_calc,
                      Paranoia  =  ExampleData[[pt]]$Paranoia[1],
                      Age       =  ExampleData[[pt]]$Age[1],
                      Sex       =  ExampleData[[pt]]$Sex[1],

                      pHI0       =  prm[1],
                      uHI0       =  prm[2],
                      pSI0       =  prm[3],
                      uSI0       =  prm[4],
                      upi        =  prm[5],
                      w0         =  prm[6],
                      wHI        =  prm[7],
                      wSI        =  prm[8],
                      eta        =  prm[9],

                      ll         =  L,
                      Success    =  conv)
          }

          fitted[[1]] <- x

        #if(simulated == 1){
          for (rep in 1:length(data)){

            evo <- infHISIll_20d(parms = as.numeric(x[rep,6:14]),
                                 d = data[[rep]][,1:3],
                                 details = T,
                                 sim_only = 1,
                                 tn = tn,
                                 phase = phase,
                                 plots = 0
                                 )

            evo <- evo$evo %>% as.data.frame() %>% dplyr::select(trial, ret, HIsim, SIsim) %>% filter(trial > 0)

            simulated[[rep]]       <- evo
            simulated[[rep]]$Order <- x[rep, 2]

          }
            simulated_df <- do.call(rbind, simulated)
            fitted[[2]] <- simulated_df
        #}

  if(plot == 1){

  fitted_plot <- ggplot(x %>%
                          dplyr::select(pHI0:eta) %>%
                          pivot_longer(1:9, 'Parameters', values_to = 'Values (Arbitary Units)'),
                        aes(`Values (Arbitary Units)`, Parameters, fill = Parameters))+

  ggridges::geom_density_ridges(stat = "binline",
                                bins = 40,
                                scale = 0.95,
                                draw_baseline = FALSE,
                                color = 'black')+
  scale_y_discrete(labels = c(
        pHI0_rec      = paste(expression(pHI[0])),
        uHI0_rec      = paste(expression(uHI[0])),
        pSI0_rec      = paste(expression(pSI[0])),
        uSI0_rec      = paste(expression(uSI[0])),
        upi_rec       = paste(expression(u[pi])),
        w0_rec        = paste(expression(w[0])),
        wHI_rec       = paste(expression(w[HI])),
        wSI_rec       = paste(expression(w[SI])),
        eta_rec       = paste(expression(eta))
  )) + theme_tq()

  likelihoodplot <- ggplot(x %>%
                             dplyr::select(ll),
                           aes(ll))+
    geom_density()+
    geom_vline(xintercept = (length(data[[1]]$Trial)*-4.38))+
    theme_tq()


  Data_Orig <- do.call(rbind, data) %>%
    as.data.frame() %>%
    plyr::join(fitted[[1]] %>%
                 as.data.frame() %>%
                 dplyr::select(ID, Order),
               by = 'ID')

  datareproduce <- ggplot(simulated_df %>%
           mutate(
               Trial = trial
           )) +

    stat_summary(data = Data_Orig, aes(Trial, HI, color = "Harmful Intent", fill = 'Harmful Intent'),  geom = "ribbon", alpha = 0.7)+
    stat_summary(data = Data_Orig, aes(Trial, HI), geom = "line", linetype = 2)+
    stat_summary(aes(Trial, HIsim), geom = "ribbon", fill  = 'grey', alpha = 0.3)+
    stat_summary(aes(Trial, HIsim), geom = "line"  )+
    stat_summary(data = Data_Orig, aes(Trial, SI, color = 'Self Interest', fill = 'Self Interest', ), geom = "ribbon", alpha = 0.7)+
    stat_summary(data = Data_Orig, aes(Trial, SI), geom = "line", linetype = 2  )+
    stat_summary(aes(Trial, SIsim), geom = "ribbon", fill  = 'grey', alpha = 0.3)+
    stat_summary(aes(Trial, SIsim), geom = "line"  )+
    labs(y = "Harmful Intent/ Self Interest")+
    scale_color_manual(values = c('#DA7000','#7700DA'), guide = 'none')+
    scale_fill_manual(values = c('#DA7000','#7700DA'), name = "Attribution")+
    facet_wrap (~ Order, scales = 'free_y', nrow = 1) +
    coord_cartesian(ylim = c(0,1))+
    tidybayes::theme_tidybayes()+
    theme(
        strip.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position = 'bottom',
        legend.direction = 'horizontal'
    )+
    theme_tq()


  plotme <- (fitted_plot | (likelihoodplot / datareproduce)) &
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.position = 'none')


    fitted[[3]] <- plotme

  }

  return(fitted)

}#end

# Utility -----------------------------------------------------------------

mysamp <- function(n, m, s, lwr, upr, nnorm) {
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}

# -   -   -   Below functions devved by Michael Moutoussis      -   -   -   -   -
# beta distribution with scaled x values, so that instead of bet. 0 and 1,
# the random var obtains values between lo and hi. Can have lo > hi, in
# which case they are mirrored. REM ncp is noncentrality parameter.
dbetasc <- function(x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE){

 xtr <- (x-lo)/(hi-lo); # will work even if hi<lo
 if (log==FALSE) {
     return( dbeta( xtr, shape1, shape2, ncp, log)/abs(hi-lo) );
 }
 else {
     return( dbeta( xtr, shape1, shape2, ncp, log) - log(abs(hi-lo)) );
 }
}
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
rbetasc <- function(n, shape1, shape2, lo=0, hi=1, ncp=0){
    auxv <- rbeta(n, shape1, shape2, ncp);
    return(  lo+ auxv*(hi-lo) );
}

# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
noisyBino <- function(pSucc, U, binN ){
   # this is the whole mdf for applying extra uncertainty U to binomial
   # distro of binN-1 draws w. success param pSucc as per
   #  MDF = binopdf(0:(binN-1),binN-1,pSuccs); MDF = MDF ^ 1/U ... etc.
 n = binN-1;
 MDF = dbinom(0:n,n,pSucc);
 MDF = MDF ^ (1/U);
 return( MDF / sum(MDF));
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
invlogit <- function(x, a=1) {
  #  inverse of logit of x bet. 0 and a
  return( a/(1+exp(-x)));
}

# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Produce k random samples from the arbirary probability vector mdf
# test with e.g. hist(rmdf(10000,dbinom(0:5,5,0.8)) )
rmdf <- function(k, mdf,check=FALSE){
    cdf <- cumsum(mdf);
    n <- length(mdf);
    if (check){
      Tol=1e-8;  # for rough checking ...
      if (abs(cdf[n] - 1) > Tol)
      { stop('mdf provided does not add up to 1 within 1e-16'); };
    }
    cdf[n] <- 1.0;  # force it to 1, for good measure ...
    x <- runif(k);
    CDF <- repRowMat(cdf,k);
    return( 1+n - rowSums(CDF >= x));
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
repRowMat <- function(x, renN) {
#
 if (!(is.vector(x))){
    stop('non - vector 1st argument in repRowMat')
  }
 return(t(matrix(rep(x,renN),nrow=length(x))));
}
