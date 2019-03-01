#High-dimensional Mediation Analysis for Survival Outcome
# load R packages
library(survival)
library(ncvreg)

# SIS for alpha
sis_alpha <- function(X, M, COV, p){
  s_alpha <- matrix(0, 3, p)
  for(j in 1:p){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    fit <- glm(M ~., data = MX)
    s_alpha[1,j] <- summary(fit)$cov.scaled[2,2]   #var for alpha
    s_alpha[2,j] <- summary(fit)$coef[2]           #coefficients for alpha
    s_alpha[3,j] <- summary(fit)$coef[2,4]         #p-value for alpha
  }
  colnames(s_alpha) = paste0("M", 1:ncol(M))
  return(s_alpha = s_alpha)
}

# SIS for beta
sis_beta <- function(X, M, Y, COV, p){
  s_beta <- matrix(0, 3, p)
  for(j in 1:p){
    if (is.null(COV)) {
      YMX <- data.frame(Y = Y, M = M[, j], X = X)
    } else {
      YMX <- data.frame(Y = Y, M = M[, j], X = X, COV = COV)
    }
    fit <- coxph(Y ~., data = YMX)
    s_beta[1,j] <- fit$var[1,1]                   #var for beta
    s_beta[2,j] <- summary(fit)$coef[1]           #coefficients for beta
    s_beta[3,j] <- summary(fit)$coef[1,5]         #p-value for beta
  }
  colnames(s_beta) = paste0("M", 1:ncol(M))
  return(s_beta = s_beta)
}

#main function
hmas <- function(X, Y, M, COV,
                 penalty = c("MCP", "SCAD", "lasso"),
                 path = c('MY', 'MX'),
                 topN = NULL,
                 verbose = TRUE, 
                 ...) {
  penalty <- match.arg(penalty)
  
  n <- nrow(M)
  p <- ncol(M)
  
  if (is.null(topN)) {
    if (path == 'MY'){
      d <- ceiling(2*n/log(n))  #the top d mediators that associated with exposure
    }else{
      d <- ceiling(3*n/log(n))
    }
  } else {
    d <- topN  
  }
  
  if(verbose) message("Step 1: Prelimenary Screening...", "     (", Sys.time(), ")")
  
  if (path == 'MY'){
    margcoef <- abs(cor(M, Y[, 1]))
    rownames(margcoef) <- paste0('M', 1:ncol(M))
    margcoef_sort <- sort(margcoef, decreasing=T)
    ID_SIS <- which(margcoef >= margcoef_sort[d])      #the index of top d mediators (Y~M)
  }else{
    alpha_s <- sis_alpha(X, M, COV, p)
    SIS_alpha <- alpha_s[2,]
    SIS_alpha_sort <- sort(SIS_alpha)
    ID_SIS <- which(SIS_alpha <= SIS_alpha_sort[d])  # the index of top d significant mediators (M~X)
  }
  
  M_SIS <- M[, ID_SIS]
  
  if(verbose) cat("Top", length(ID_SIS), "mediators selected (ranked from most to least significant): ", colnames(M_SIS), "\n")
  
  XM <- cbind(M_SIS, X)
  C_M<-colnames(XM)
  
  
  if(verbose) message("Step 2: Penalized Variable Selection (", penalty, ") ...", "   s  (", 
                      Sys.time(), ")")
  
  if (is.null(COV)) {
    fit <- ncvsurv(XM, Y, 
                   penalty = penalty, 
                   penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...)
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1,drop=F]
    conf.names <- colnames(COV)
    XM_COV <- cbind(XM, COV)
    fit <- ncvsurv(XM_COV, Y, 
                   penalty = penalty, 
                   penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV))), ...)
  }
  
  lam <- fit$lambda[which.min(BIC(fit))]
  if(verbose) cat("lambda selected: ", lam, "\n")
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[1:d]
  ID_p_non <- which(est != 0)
  
  if(verbose) cat("Non-zero", penalty, "beta estimate(s) of mediator(s) found: ", names(ID_p_non), "\n")
  
  beta_p <- est[ID_p_non]  # The non-zero MCP estimators of beta
  ID_p <- ID_SIS[ID_p_non]  # The index of the ID of non-zero beta in the Cox regression
  MCP_M <- names(ID_p_non)
  
  
  if(verbose) message("Step 3: The adjusted Sobel significance test ...", 
                      "     (", Sys.time(), ")")
  
  if (path == 'MY'){
    sis_alpha <- sis_alpha(X, M, COV, p)
    alpha_s <- sis_alpha[, ID_p]
  }else{
    alpha_s <- alpha_s[, ID_p]
  }
  
  alpha_est <- alpha_s[2, ]   #  the estimator for alpha
  var_alpha <- alpha_s[1, ]
  
  if (is.null(COV)) {
    YMX <- data.frame(Y = Y, M[, ID_p, drop = FALSE], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[, ID_p, drop = FALSE], X = X, COV = COV)
  }
  
  cox_model <- coxph(Y ~ ., data = YMX)
  
  beta_est <- summary(cox_model)$coefficients[1: length(ID_p)]     #the estimator of beta
  DE <- summary(cox_model)$coefficients[(length(ID_p)+1), 1]
  DE <- exp(DE)
  
  var_beta <- diag(cox_model$var[1:length(ID_p),1:length(ID_p)])
  
  ab_est <- alpha_est * beta_est   # the estimator of alpha*beta
  
  # var(alpha*beta)
  var_ab <- (alpha_est^2) * var_beta + var_alpha * (beta_est^2)
  
  # confidence interval
  conf_low <- ab_est - 1.96 * sqrt(var_ab)
  conf_up <- ab_est + 1.96 * sqrt(var_ab)
  
  # sobel test for alpha*beta
  s.test <- abs(ab_est)/sqrt(var_ab)   #z-score for sobel test
  P_sobel <- 2 * (1-pnorm(s.test))     #p-value of sobel test
  P_bon_sobel <- p.adjust(P_sobel, 'bonferroni', length(ID_p)) #Bonferroni adjusted p-value
  P_bon_sobel[P_bon_sobel > 1] <- 1
  P_fdr_sobel <- p.adjust(P_sobel, 'fdr', length(ID_p)) #FDR adjusted p-value
  P_fdr_sobel[P_fdr_sobel > 1] <- 1
  Pf_sobel <- as.matrix(P_fdr_sobel)
  rownames(Pf_sobel) <- paste0('M',ID_p)
  ID.a <- rownames(Pf_sobel)
  ID_test <- which(Pf_sobel < 0.05)
  ID_t <- ID.a[ID_test]
  
  if(verbose) cat("Significant", "mediator(s) found: ", ID_t, "\n")
  
  results <- data.frame(alpha = alpha_est, beta = beta_est,
                        `alpha_est*beta_est` = ab_est, 
                        conf_low=conf_low, conf_up=conf_up, s.test=s.test, 
                        P_sobel=P_sobel, P_bon_sobel=P_bon_sobel, P_fdr_sobel=P_fdr_sobel,
                        var_ab=var_ab, var_alpha=var_alpha, var_beta=var_beta,
                        DE, check.names = FALSE)
  
  
  if(verbose) message("Done!", "     (", Sys.time(), ")")
  
  return(list(C_M, MCP_M, ID_t, results))
}

#run the model
hmas.fit_M <- hmas(X=X, Y=Y, M=M, COV=COV, path='MX', verbose=TRUE) 
