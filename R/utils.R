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