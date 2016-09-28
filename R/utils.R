test2groups <- function(df, varlist, groupvarname = "group", group1 = "P", group2 = "m", meansd_decpoints = 2, pval_decpoints = 3) {
  n <- length(varlist)
  dfout <- data.frame(variable = varlist, group1mean = numeric(n), group1sd = numeric(n), group2mean = numeric(n), group2sd = numeric(n), t = numeric(n), df = numeric(n), p = numeric(n), d = numeric(n), d_CI = character(n), stringsAsFactors = F)
  for (i in 1:n) {
    g2 <- eval(parse(text = paste0("df[df$",groupvarname, "=='",group1,"','", varlist[i],"']")))[[1]]
    g1 <- eval(parse(text = paste0("df[df$",groupvarname, "=='",group2,"','", varlist[i],"']")))[[1]]

    dfout$group1mean[i] <- g1 %>% mean(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group1sd[i] <- g1 %>% sd(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group2mean[i] <- g2 %>% mean(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group2sd[i] <- g2 %>% sd(na.rm = T) %>% round(digits = meansd_decpoints)
    tt <- t.test(g1, g2, var.equal = T)

    dfout$t[i] <- round(tt$statistic, meansd_decpoints)
    dfout$df[i] <- tt$parameter
    dfout$p[i] <- round(tt$p.value, pval_decpoints)

    dd <- compute.es::tes(tt$statistic, length(g2), length(g1), verbose = F)

    dfout$d[i] <- round(dd$d, meansd_decpoints)
    dfout$d_CI[i] <- paste0("[", dd$l.d, ", ", dd$u.d, "]")

  }
  colnames(dfout)[2:5] <- c(sprintf("group%smean", group2), sprintf("group%ssd", group2),
                            sprintf("group%smean", group1), sprintf("group%ssd", group1))
  dfout
}

test3groups <- function(df, varlist, groupvarname = "group", group1, group2, group3, meansd_decpoints = 2, pval_decpoints = 3, boot_CI = F) {
  n <- length(varlist)
  dfout <- data_frame(variable = varlist, group1mean = numeric(n), group1sd = numeric(n), group2mean = numeric(n), group2sd = numeric(n), group3mean = numeric(n), group3sd = numeric(n), F = numeric(n), df1 = numeric(n), df2 = numeric(n), p = numeric(n), eta_sq = numeric(n), eta_sq_CI = character(n))
  
  for (i in 1:n) {
    g1 <- eval(parse(text = paste0("df[df$",groupvarname, "=='",group1,"','", varlist[i],"']")))[[1]]
    g2 <- eval(parse(text = paste0("df[df$",groupvarname, "=='",group2,"','", varlist[i],"']")))[[1]]
    g3 <- eval(parse(text = paste0("df[df$",groupvarname, "=='",group3,"','", varlist[i],"']")))[[1]]
    
    dfout$group1mean[i] <- g1 %>% mean(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group1sd[i] <- g1 %>% sd(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group2mean[i] <- g2 %>% mean(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group2sd[i] <- g2 %>% sd(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group3mean[i] <- g3 %>% mean(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group3sd[i] <- g3 %>% sd(na.rm = T) %>% round(digits = meansd_decpoints)
    
    dfG <- data.frame(val = c(g1,g2,g3),group = c(rep(group1,length(g1)), rep(group2,length(g2)), rep(group3,length(g3))))
    dfG$id <- 1:nrow(dfG)
    aov1 <- ez::ezANOVA(dfG[complete.cases(dfG),],dv = .(val), between = .(group), wid = .(id))
    
    dfout$F[i] <- round(aov1$ANOVA$F, meansd_decpoints)
    dfout$df1[i] <- aov1$ANOVA$DFn
    dfout$df2[i] <- aov1$ANOVA$DFd
    dfout$p[i] <- round(aov1$ANOVA$p, pval_decpoints)
    
    if (boot_CI) {
      aov_etasq_p <- anova_etasq_p_ci(dfG[complete.cases(dfG),], formula(val ~ group))
    } else {
    dfout$eta_sq[i] <- round(aov1$ANOVA$ges, meansd_decpoints)
    
    ci <- MBESS::ci.pvaf(F.value = aov1$ANOVA$F,
            df.1 = dfout$df1[i],
            df.2 = dfout$df2[i],
            N = nrow(dfG[complete.cases(dfG),]),
            conf.level = .95)
    
    dfout$eta_sq_CI[i] <- paste0("[", ci$Lower.Limit.Proportion.of.Variance.Accounted.for %>% round(meansd_decpoints), ", ", ci$Upper.Limit.Proportion.of.Variance.Accounted.for %>% round(meansd_decpoints), "]")
    
    }
  }
  colnames(dfout)[2:7] <- c(sprintf("group%smean", group1), sprintf("group%ssd", group1),
                            sprintf("group%smean", group2), sprintf("group%ssd", group2),
                            sprintf("group%smean", group3), sprintf("group%ssd", group3))
  dfout
}

etasq_p <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  dfx <- lm(formula, data = d) %>% car::Anova() %>% broom::tidy()
  
  ressq <- dfx %>% filter(term == "Residuals") %>% .$sumsq
  eta.sq <- dfx %>% group_by(term) %>% mutate(eta.sq = sumsq / (sumsq + ressq)) %>% .$eta.sq
  return(eta.sq)
} 


anova_etasq_p_ci <- function(df, f, nSamp = 1000) {
  dfx <- lm(f, data = df) %>% car::Anova() %>% broom::tidy()
  dfx$eta.sq.p <- NA
  dfx$eta.sq.p.lo <- NA
  dfx$eta.sq.p.hi <- NA
  
  results <- boot::boot(data = df, statistic = etasq_p, R=nSamp, formula = f)
  
  for (i in 1:(length(results$t0)-1)){
    dfx$eta.sq.p[i] <- results$t0[i]
    
    bci <- boot::boot.ci(results, type="bca",index = i) 
    dfx$eta.sq.p.lo[i] <- bci$bca[length(bci$bca) - 1]
    dfx$eta.sq.p.hi[i] <- bci$bca[length(bci$bca)]
  }
  dfx <- dfx %>% mutate(statistic = round(statistic,2),p.value = round(p.value,3),eta.sq.p = round(eta.sq.p,2), eta.sq.p.lo = round(eta.sq.p.lo,2), eta.sq.p.hi = round(eta.sq.p.hi,2))
  return(dfx)
  
  
}

#' normalize_flow
#'
#' @param df - data frame with flow scores
#' @param flowcol - name of the column containg flow scores
#'
#' @return data frame with flow scores standardized by norms  
#' @export
#'
#' @examples
normalize_flow <- function(df, flowcol) {
  df2 <- df %>% left_join(flow_norms, by = setNames("HS",flowcol))
  df2[[flowcol]] <- df2$tskor
  df <- df2 %>% select(-tskor)
  return(df)
}

#' Performs many one sample tests and outputs as data_frame()
#'
#' @param df - data frame containing given variables
#' @param varlist - variables to be tested
#' @param testvals - mean values to be tested against. If only one value is used, all variables will be tested against this value
#'
#' @return data frame with results for one sample test
#' @export
#'
#' @examples
multiple_onesample_ttests <- function(df, varlist, testvals, meansd_decpoints = 2, pval_decpoints = 3) {
  if (length(testvals) == 1) {
    testvals <- rep(testvals, times = length(varlist))
  }
  
  n <- length(varlist)
  
  stopifnot(length(testvals) == n)
  
  dfout <- data_frame(variable = varlist, tested_val = testvals, mean = numeric(n), sd = numeric(n), t = numeric(n), df = numeric(n), p = numeric(n))
  
  for (i in 1:n) {
    g1 <- df[[varlist[i]]]
    
    
    dfout$mean[i] <- g1 %>% mean(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$sd[i]   <- g1 %>% sd(na.rm = T) %>% round(digits = meansd_decpoints)
    tt <- t.test(g1, mu = testvals[i])
    
    dfout$t[i] <- round(tt$statistic, meansd_decpoints)
    dfout$df[i] <- tt$parameter
    dfout$p[i] <- round(tt$p.value, pval_decpoints)
  }
  return(dfout)
}
