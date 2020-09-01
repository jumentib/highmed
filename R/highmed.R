##' Epigenome Wide Association Study with both exposure and outcome
##'
##' This function uses lfmm (latent factor mixed models) to estimate
##' the effects of exposures and outcomes on a response matrix.
##'
##'
##' @param M a response variable matrix with n rows and p columns.
##' Each column corresponds to a beta-normalized methylation profile.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param K an integer for the number of latent factors in the regression model.
##' @param covar set of covariable, must be numeric.
##' @return an object with the following attributes:
##'
##'  - U the latent variable score matrix with dimensions n x K.
##'
##'  - B the effect size matrix for the exposure X and the outcome Y.
##'
##'  - score matrix for the exposure X and the outcome Y.
##'
##'  - pValue matrix for the exposure X and the outcome Y.
##'
##'  - calibrated.score2, the calibrated score matrix for the exposure X and the outcome Y.
##'
##'  - calibrated.pvalue, the calibrated pValue matrix for the exposure X and the outcome Y.
##'
##'  - GIF : Genomic Inflation Factor for exposure and outcome
##'
##' @details
##' The response variable matrix Y and the explanatory variable are centered.
##' Missing values must be imputed. The number of latent factors can be estimated
##' by looking at the screeplot of eigenvalues of a PCA.
##' Possibility of calibrating the scores and pValues by the GIF (Genomic Inflation Factor).
##' See lfmm package for more information.
##' @export
##' @author Basile Jumentier
##' @references
##' @examples
##'
multivariate_EWAS <- function(X, Y, M, K, covar = NULL) {

  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, Y, covar), K = K)
  beta <- dat$B[, 1:2]
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, Y, covar), lfmm = dat)

  return(list(beta = beta,
              score = dat$score[, 1:2],
              pValue = dat$pvalue[, 1:2],
              calibrated.score2 = dat$calibrated.score2[, 1:2],
              calibrated.pvalue = dat$calibrated.pvalue[, 1:2],
              gif = dat$gif))
}



##' Compute the squared maximum of two series of pValues
##'
##' This function compute the squared maximum of two series of pValues from the multivariate_EWAS() function.
##' The objective of this function is to test all the markers and to determine which could be
##' potential mediators in the exposure-outcome association.
##'
##' @param pval1 vector of pValues (p*1) of exposure.
##' @param pval2 vector of pValues (p*1) of ouctome.
##' @param diagnostic.plot if TRUE the histogram of the p-values together
##' with the estimate of eta0 null line is plotted.
##' This is useful to visually check the fit of the estimated proportion of null p-values.
##' @param ... argument of the fdrtool function from the fdrtool package
##'
##' @return an object with the following attributes:
##'
##'  - a pValue for each markers
##'
##'  - a qValue for each markers
##'
##'  - the eta0 of the set of pValues
##'
##' @details
##' The pValue is computed for each markers following this formula
##'
##' \deqn{pV = max(pVal1, pVal2)^2}
##'
##' This quantity eta0, i.e. the proportion eta0 of null p-values in a given vector of p-values,
##' is an important parameter when controlling the false discovery rate (FDR).
##' A conservative choice is eta0 = 1 but a choice closer to the true value will
##' increase efficiency and power - see Benjamini and Hochberg (1995, 2000) and Storey (2002) for details.
##' We use the fdrtool package to transform pValues into qValues,
##' which allows us to control the FDR.
##' @export
##' @author Basile Jumentier
##' @references
##' @examples
##'
max2 <- function(pval1, pval2, diagnostic.plot = F, ...) {

  pval <- apply(cbind(pval1, pval2), 1, max)^2
  eta0 <- fdrtool::pval.estimate.eta0(pval, diagnostic.plot = diagnostic.plot)
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", plot = F, verbose = F, ...)$qval

  return(list(pval = pval,
              eta0 = eta0,
              qval = qval))
}

##' Run mediation analysis for a set of markers
##'
##' Estimate various quantities for causal mediation analysis for each
##' significant markers, including average causal mediation effects
##' (indirect effect), average direct effects, proportions mediated,
##' and total effect.
##'
##' @param qval set of qValues from max2() function
##' @param M a response variable matrix with n rows and p columns.
##' Each column corresponds to a beta-normalized methylation profile.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param covar set of covariable, must be numeric.
##' @param U set of latent factors from multivariate_EWAS() function
##' @param FDR FDR threshold to pass markers in mediation analysis
##' @param sims number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
##' 10000 is recommended.
##' @param ... argument of the mediate function from the mediation package
##'
##' @return
##' Tables of results of mediation analyzes for markers with a qValue below the FDR threshold.
##' Indirect effect (ACME - average causal mediation effect), ADE (average direct effect),
##' PM (proportion mediated) and TE (total effect). Composition of tables: estimated effect,
##' confidence interval and mediation pValue.
##'
##' @details
##'
##' We use the mediate() function of the mediation package on the set of markers having a qValue lower
##' than the FDR threshold. This function makes it possible to estimate their indirect effects and to
##' test their significance.
##'
##' @export
##' @author Basile Jumentier
##' @references
##' @examples
##'
univariate_mediation <- function(qval, X, Y, M, covar = NULL, U = NULL, FDR = 0.1, sims = 3, ...) {

  if (is.null(colnames(M))) {
    colnames(M) <- 1:ncol(M)
  }

  M <- M[, qval <= FDR]


  ACME <- matrix(ncol = 4, nrow = ncol(M))
  ADE <- matrix(ncol = 4, nrow = ncol(M))
  PM <- matrix(ncol = 4, nrow = ncol(M))
  TE <- matrix(ncol = 4, nrow = ncol(M))

  for (i in 1:ncol(M)) {

    dat.x <- data.frame(X = X, Mi = M[, i], covar = cbind(covar, U))
    dat.y <- data.frame(X = X, Mi = M[, i], covar = cbind(covar, U), Y = Y)

    mod1 <- lm(Mi ~ X + ., data = dat.x)
    mod2 <- lm(Y ~ X + Mi + ., data = dat.y)

    med <- mediation::mediate(mod1, mod2, sims = sims, treat = "X", mediator = "Mi", ...)

    ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
    ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
    PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
    TE[i, ] <- c(med$tau.coef, med$tau0.ci[1], med$tau0.ci[2], med$tau.p)
  }

  ACME <- as.data.frame(ACME)
  ADE <- as.data.frame(ADE)
  PM <- as.data.frame(PM)
  TE <- as.data.frame(TE)

  colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")

  ACME$CpG <- colnames(M)
  ADE$CpG <- colnames(M)
  PM$CpG <- colnames(M)
  TE$CpG <- colnames(M)

  return(list(ACME = ACME,
              ADE = ADE,
              PM = PM,
              TE = TE))

}

##' Summary plot for ACME
##'
##' This function draw a summary plot of ACME (average causal mediation effect)
##'
##' @param ACME the table of ACME from the univariate_mediation() function
##' @return
##' Summary plot for ACME
##'
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
plot_summary_ACME <- function(ACME) {

  p <- ggplot(ACME, aes(est, reorder(CpG, est), color = pval <= 0.05, shape = pval <= 0.05)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5)) +
    geom_point(size = 0.8) +
    theme_bw() +
    xlab("ACME (Average Causal Mediation Effect)") +
    ylab("CpG") +
    theme(panel.border = element_blank(),
          panel.spacing = unit(0.01, "lines"),
          axis.ticks = element_blank()) +
    scale_color_manual(values = c("black", "red"))

  print(p)
}



##' Summary plot for univariate_mediation function
##'
##' This function draw a summary plot of the mediation analysis
##'
##' @param res_univariate_mediation result object from univariate_mediation() function
##' @return
##' Summary plot
##'
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
plot_summary_med <- function(res_univariate_mediation) {

  tmp <- rbind(cbind(res$ACME, stat = "ACME"),
               cbind(res$ADE, stat = "ADE"),
               cbind(res$PM, stat = "PM"),
               cbind(res$TE, stat = "TE"))

  p <- ggplot(tmp, aes(est, stat, color = stat, shape = stat)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5)) +
    geom_point() +
    facet_grid(CpG ~ ., scales = "free_y") +
    theme_bw() +
    xlab(NULL) +
    ylab(NULL) +
    labs(color = NULL, shape = NULL) +
    theme(strip.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "grey"),
          strip.text.y = element_text(angle = 0),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.x = unit(0.1, "lines"),
          panel.spacing.y = unit(0.01, "lines"),
          legend.position = "left") +
    scale_color_manual(values = c("red", "blue", "orange", "black"))

  print(p)
}

##' Run mediation analysis for a set of markers
##'
##' Function adapt from the combp function() of the ENmix package
##'
##' @param data A data frame from bed format file with colname name
##' "V1","V2", "V3","V4","V5",V1 indicate chromosome (1,2,3,...,X,Y),
##' V2 is chromosome position, V4 is for P value and V5 for name of CpGs.
##' @param dist.cutoff Maximum distance in base pair to combine adjacent DMRs.
##' @param bin.size bin size for autocorrelation calculation.
##' @param seed FDR significance threshold for initial selection of DMR region.
##' @param nCores Number of computer cores used in calculation
##'
##' @return
##' Results of the DMRs analysis.
##'
##' @details
##'
##' The input should be a data frame with column name V1-V5, indicating chromosome, start position,end position,
##' pValues and probe names. The function will use a modified comb-p method to identify
##' differentially methylated regions.
##'
##' @author Basile Jumentier
##' @references
##' @examples
##'
combp2 <- function (data, dist.cutoff = 1000, bin.size = 310, seed = 0.01, nCores = 10) {

  ##### a function to get a table of p-values for estimating acf
  #####loc should be increasing;
  acf.table<-function(x,loc,dist.cutoff){
    flag=TRUE; lag=1; result=NULL
    while(flag){
      x1=utils::head(x,-lag); x2=utils::tail(x,-lag); dist=diff(loc,lag=lag)
      index=(dist<dist.cutoff)
      if(all(!index)){flag=FALSE}else{
        result=rbind(result,data.frame(x1=x1[index],x2=x2[index],dist=dist[index]))
        lag=lag+1
      }
    }
    return(result)
  }

  ##### a function to estimate acf
  get.acf<-function(data,dist.cutoff,bin.size){
    temp<-NULL
    for (chr in unique(data$V1)){
      y<-data[data$V1==chr,]; y<-y[order(y$V3),]
      temp<-rbind(temp,acf.table(y$V4,y$V3,dist.cutoff))
    }
    bin.label<-findInterval(temp$dist,seq(bin.size,dist.cutoff,bin.size))
    temp.stouffer<-by(temp,bin.label,FUN=function(x){cor.test(qnorm(x$x1),
                                                              qnorm(x$x2),alternative="greater")},simplify=FALSE)

    cor.stouffer<-sapply(temp.stouffer,function(x){x$estimate})
    p.stouffer<-sapply(temp.stouffer,function(x){x$p.value})

    if (any(p.stouffer>0.05)){
      index=min(which(p.stouffer>0.05))
      cor.stouffer[index:length(cor.stouffer)]=0
    }
    return(cor.stouffer)
  }

  if (nCores > detectCores()) {
    nCores = detectCores()
  }
  data = as.data.frame(data)
  acf <- get.acf(data, dist.cutoff, bin.size)
  result <- mclapply(unique(data$V1), function(chr) {
    y = data[data$V1 == chr, ]
    y = y[order(y$V3), ]
    pos = y$V3
    p = qnorm(y$V4)
    temp = sapply(pos, function(i) {
      index.i = (abs(pos - i) < bin.size)
      if (sum(index.i) > 1) {
        int <- findInterval(c(dist(pos[index.i])), c(bin.size,
                                                     2 * bin.size))
        sd <- sqrt(sum(acf[int + 1]) * 2 + sum(index.i))
        return(pnorm(sum(p[index.i]), mean = 0, sd = sd))
      }
      else {
        return(y$V4[index.i])
      }
    })
    return(data.frame(chr, start = pos, end = pos, s.p = temp))
  }, mc.cores = nCores)
  result <- do.call("rbind", result)
  names(result) = c("chr", "start", "end", "s.p")
  result = result[p.adjust(result$s.p, method = "fdr") < seed,
                  ]
  result.fdr = NULL
  if (nrow(result) > 0) {
    for (chr in unique(result$chr)) {
      y = data[data$V1 == chr, ]
      y = y[order(y$V3), ]
      pos = y$V3
      p = qnorm(y$V4)
      result.chr = result[result$chr == chr, ]
      a = IRanges::IRanges(start = result.chr$start, end = result.chr$end)
      b = IRanges::reduce(a, min.gapwidth = dist.cutoff)
      start = start(b)
      end = end(b)
      region.max <- max(width(b))
      temp = sapply(1:length(b), function(i) {
        index.i = (pos >= start[i] & pos <= end[i])
        if (sum(index.i) > 1) {
          int <- findInterval(c(dist(pos[index.i])),
                              seq(bin.size, region.max + bin.size, bin.size))
          sd <- sqrt(sum(ifelse(int < length(acf), acf[int +
                                                         1], 0)) * 2 + sum(index.i))
          return(pnorm(sum(p[index.i]), mean = 0, sd = sd))
        }
        else {
          return(y$V4[index.i])
        }
      })
      result.fdr = rbind(result.fdr, data.frame(chr, start,
                                                end, p = temp))
    }
    result.fdr$fdr = p.adjust(result.fdr$p, method = "fdr")
    result.fdr <- result.fdr[order(result.fdr$p), ]
    result.fdr$start = (result.fdr$start - 1)
  }

  return(result.fdr)
}


##' Find DMR
##'
##' To identify differentially methylated regions using a modified comb-p method.
##'
##' @param chr chromosomes
##' @param start chromosomal position of markers (start)
##' @param end chromosomal position of markers (end)
##' @param pval pValues for each markers, from the max2 function
##' @param cpg name of each markers
##' @param ... argument of the combp function from the ENmix package
##'
##' @return
##' A set of selected DMRs.
##'
##' @details
##'
##' The function will use a modified comb-p method to identify
##' differentially methylated regions (DMRs).
##'
##' @export
##' @author Basile Jumentier
##' @references
##' @examples
##'
DMR_search <- function(chr, start, end, pval, cpg, ...) {

  tmp <- data.frame(chr, start, end, pval, cpg)
  colnames(tmp) <- paste0("V", 1:5)

  tmp <- combp2(tmp, ...)

  return(list(res = tmp,
              data = data.frame(chr, start, end, pval, cpg)))
}


##' Build DMR vector
##'
##' To build a vector for each DMR find with DMR_search
##'
##' @param res result object of DMR_search function
##' @param methylation chromosomal position of markers (start)
##' @param nb_cpg chromosomal position of markers (end)
##'
##' @return
##' A set of build DMRs.
##'
##' DMR_acp contains the first components of each PCA for each DMR.
##' CpG_for_each_DMR contains the list of markers (CpGs) present on each DMR.
##'
##'
##' @details
##'
##' We use the series of pValues (one pValue per CpGs) obtained with the EWASmultivariate
##' regression method and the combination of pValue max2.
##' To determine the potential DMRs used the combp method present in the ENmix package (Xu et al. 2016).
##' This method uses the Fisher method to combine the pValues and also the base pair distance (bP)
##' between CpGs (1000 bP maximum between nb_cpg CpGs on the same DMR).
##' The information for each DMR is summarized by running a PCA by DMR on all of the CpGs present on each DMR.
##' Recovering the first principal component of each PCA,
##' we therefore have a vector corresponding to the first principal component of a PCA for each DMR.
##'
##' @export
##' @author Basile Jumentier
##' @references
##' @examples
##'
DMR_built <- function(res, methylation, nb_cpg = 2) {

  data <- res$data
  res <- res$res

  # Number of CpG per DMR

  nb <- NULL

  for (i in 1:nrow(res)) {

    chri <- as.character(res$chr[i])

    tmp <- dplyr::filter(data, chr == chri)

    nb <- c(nb, sum((res$start[i]:res$end[i]) %in% tmp$start))
  }

  # Select DMRs with nb_cpg CpGs at minimum

  res <- cbind(res, nb)

  res <- dplyr::filter(res, nb >= nb_cpg)

  DMR.select <- list()

  for (i in 1:nrow(res)) {

    chri <- as.character(res$chr[i])

    tmp <- dplyr::filter(data, chr == chri)

    DMR.select[[i]] <- tmp$cpg[(tmp$start %in% (res$start[i]:res$end[i]))]
  }

  # Select CpGs values in the methylation matrix

  DMR.meth <- list()

  for (i in 1:length(DMR.select)) {
    DMR.meth[[i]] <- methylation[, DMR.select[[i]]]
  }

  # Built a vector for each DMR with the first component of PCA

  DMR.acp <- as.data.frame(matrix(ncol = length(DMR.meth), nrow = nrow(methylation)))
  colnames(DMR.acp) <- paste0("DMR", 1:length(DMR.meth))

  for (i in 1:length(DMR.meth)) {
    DMR.acp[, i] <- prcomp(DMR.meth[[i]])$x[, 1]
  }

  # data

  res <- cbind(DMR = colnames(DMR.acp), res)
  names(DMR.select) <- colnames(DMR.acp)

  return(list(DMR_acp = DMR.acp,
              res = res,
              CpG_for_each_DMR = DMR.select))
}

##' Run mediation analysis on DMR
##'
##' Estimate various quantities for causal mediation analysis for each
##' DMRs, including average causal mediation effects
##' (indirect effect), average direct effects, proportions mediated,
##' and total effect.
##'
##' @param DMR a matrix of DMRs from the DMR_built() function (DMR_acp).
##' @param X an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param covar set of covariable, must be numeric.
##' @param U set of latent factors from multivariate_EWAS() function
##' @param sims number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
##' 10000 is recommended.
##' @param ... argument of the mediate function from the mediation package
##'
##' @return
##' Tables of results of mediation analyzes for markers with a qValue below the FDR threshold.
##' Indirect effect (ACME - average causal mediation effect), ADE (average direct effect),
##' PM (proportion mediated) and TE (total effect). Composition of tables: estimated effect,
##' confidence interval and mediation pValue.
##'
##' @details
##'
##' We use the mediate() function of the mediation package on the set of selected DMRs.
##' This function makes it possible to estimate their indirect effects and to
##' test their significance.
##'
##' @export
##' @author Basile Jumentier
##' @references
##' @examples
##'
univariate_mediation_DMR <- function(X, Y, DMR, covar = NULL, U = NULL, sims = 3) {

  ACME <- matrix(ncol = 4, nrow = ncol(DMR))
  ADE <- matrix(ncol = 4, nrow = ncol(DMR))
  PM <- matrix(ncol = 4, nrow = ncol(DMR))
  TE <- matrix(ncol = 4, nrow = ncol(DMR))

  for (i in 1:ncol(DMR)) {

    dat.x <- data.frame(X = X, Mi = DMR[, i], covar = cbind(covar, U))
    dat.y <- data.frame(X = X, Mi = DMR[, i], covar = cbind(covar, U), Y = Y)

    mod1 <- lm(Mi ~ X + ., data = dat.x)
    mod2 <- lm(Y ~ X + Mi + ., data = dat.y)

    med <- mediation::mediate(mod1, mod2, sims = sims, treat = "X", mediator = "Mi")

    ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
    ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
    PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
    TE[i, ] <- c(med$tau.coef, med$tau0.ci[1], med$tau0.ci[2], med$tau.p)
  }

  ACME <- as.data.frame(ACME)
  ADE <- as.data.frame(ADE)
  PM <- as.data.frame(PM)
  TE <- as.data.frame(TE)

  colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")

  ACME$DMR <- colnames(DMR)
  ADE$DMR <- colnames(DMR)
  PM$DMR <- colnames(DMR)
  TE$DMR <- colnames(DMR)

  return(list(ACME = ACME,
              ADE = ADE,
              PM = PM,
              TE = TE))

}




##' Summary plot for univariate_mediation_DMR function
##'
##' This function draw a summary plot of the mediation analysis
##'
##' @param res_univariate_mediation_DMR result object from res_univariate_mediation_DMR() function
##' @param res_DMR_built result object from res_DMR_built() function
##' @return
##' Summary plot
##'
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
plot_summary_DMR <- function(res_univariate_mediation_DMR, res_DMR_built) {

  tmp <- merge.data.frame(res_univariate_mediation_DMR$ACME,
                          res_DMR_built$res, by.x = 5, by.y = 1)

  tmp$dmr <- paste0(tmp$chr, ":", tmp$start, "-", tmp$end)

  p <- ggplot(tmp, aes(est, reorder(dmr, pval), color = nb, shape = pval <= 0.05)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point() +
    geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5)) +
    theme_bw() +
    xlab("ACME (Average Causal Mediation Effect)") +
    ylab("DMRs") +
    labs(color = "Number of \nCpGs on DMRs",
         shape = "pValues of \nmediation <= 0.05") +
    theme(panel.border = element_blank(),
          panel.spacing = unit(0.01, "lines"),
          axis.ticks = element_blank()) +
    scale_color_gradient(low = "blue", high = "red")


  print(p)
}
