library(data.table); library(logistf)
co  <- readRDS("out/02_cohort_GSE206285.rds")
ex  <- co$expr_all
arm <- factor(co$arm, levels = c("placebo", "treated"))
y   <- co$y
g   <- intersect(c("OSM", "OSMR"), rownames(ex)); print(g)
z   <- t(scale(t(ex[g, , drop = FALSE])))
sc  <- setNames(lapply(g, function(x) z[x, ]), g)
if (length(g) == 2) sc$OSM_axis <- colMeans(z)

res <- rbindlist(lapply(names(sc), function(n) {
  d <- data.frame(y = y, arm = arm, sc = as.numeric(scale(sc[[n]])))
  m <- logistf(y ~ arm * sc, data = d)
  k <- grep(":", names(coef(m))); j <- which(names(coef(m)) == "sc")
  data.table(marker = n,
             OR_placebo_arm = exp(coef(m)[j]), lo = exp(m$ci.lower[j]),
             hi = exp(m$ci.upper[j]), p_placebo = m$prob[j],
             inter_OR = exp(coef(m)[k]), inter_lo = exp(m$ci.lower[k]),
             inter_hi = exp(m$ci.upper[k]), inter_p = m$prob[k])
}))
print(as.data.frame(res), digits = 3)
fwrite(res, "out/06c_west_OSM.csv")
