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

test3groups <- function(df, varlist, groupvarname = "group", group1, group2, group3, meansd_decpoints = 2, pval_decpoints = 3) {
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
    dfout$group3mean[i] <- g2 %>% mean(na.rm = T) %>% round(digits = meansd_decpoints)
    dfout$group3sd[i] <- g2 %>% sd(na.rm = T) %>% round(digits = meansd_decpoints)
    
    dfG <- data.frame(val = c(g1,g2,g3),group = c(rep(group1,length(g1)), rep(group2,length(g2)), rep(group3,length(g3))))
    dfG$id <- 1:nrow(dfG)
    aov1 <- ez::ezANOVA(dfG[complete.cases(dfG),],dv = .(val), between = .(group), wid = .(id))
    
    dfout$F[i] <- round(aov1$ANOVA$F, meansd_decpoints)
    dfout$df1[i] <- aov1$ANOVA$DFn
    dfout$df2[i] <- aov1$ANOVA$DFd
    dfout$p[i] <- round(aov1$ANOVA$p, pval_decpoints)
    
    dfout$eta_sq[i] <- round(aov1$ANOVA$ges, meansd_decpoints)
    
    ci <- MBESS::ci.pvaf(F.value = dfout$F[i],
            df.1 = dfout$df1[i],
            df.2 = dfout$df2[i],
            N = nrow(dfG[complete.cases(dfG),]),
            conf.level = .95)
    
    dfout$eta_sq_CI[i] <- paste0("[", ci$Lower.Limit.Proportion.of.Variance.Accounted.for %>% round(meansd_decpoints), ", ", ci$Upper.Limit.Proportion.of.Variance.Accounted.for %>% round(meansd_decpoints), "]")
    
    
  }
  colnames(dfout)[2:7] <- c(sprintf("group%smean", group1), sprintf("group%ssd", group1),
                            sprintf("group%smean", group2), sprintf("group%ssd", group2),
                            sprintf("group%smean", group3), sprintf("group%ssd", group3))
  dfout
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
