# Extract the additional Table 2 rows from the same estimators used for the
# published analysis; import the forest data preparation, never its plotting.
root <- normalizePath(".")
.libPaths(c(file.path(root,"revision/R_lib"),.libPaths()))
setwd(file.path(root,"meta_analysis"))
s <- readLines("01_main_analysis.R",encoding="UTF-8")
s[grepl("^setwd\\(",s)] <- ""
invisible(capture.output(eval(parse(text=s),envir=.GlobalEnv)))
out <- list()
add <- function(label,d,fit) {
  out[[length(out)+1L]] <<- data.frame(analysis=label,effects=nrow(d),
    clusters=length(unique(d$dep_id)),studies=length(unique(d$study_id)),
    estimate=fit$pooled_kappa,ci_lower=fit$ci_kappa[1],ci_upper=fit$ci_kappa[2],
    pi_lower=fit$pi_kappa[1],pi_upper=fit$pi_kappa[2])
}
for (x in list(c("ir_wk","IR weighted kappa (+AT-derived LW)"),c("im_wk","IM weighted kappa (+AT-derived LW)"),
  c("ir_uk","IR unweighted kappa"),c("im_uk","IM unweighted kappa"),c("ir_5c","IR 5-cat kappa"),c("im_5c","IM 5-cat kappa"))) add(x[2],get(x[1]),get(paste0(x[1],"_res")))
add("Intra-reader kappa all-effects",intra_all_k,intra_rve)
add("Inter-version kappa all-effects",iv_all_k,iv_rve)
for (q in parse("04b_forest_compact.R",encoding="UTF-8"))
  if(is.call(q) && identical(q[[1]],as.name("<-")) && identical(q[[2]],as.name("prep_opes"))) eval(q)
op <- prep_opes(comp %>% filter(is.na(exclude)|exclude=="",opes_include==1,!is.na(as.numeric(kappa))))
ir <- subset(op,analytic_stratum=="inter-reader")
im <- subset(op,analytic_stratum=="inter-modality")
groups <- list(IR=ir,IM=im,CT=subset(ir,mod_group=="CT"),MRI=subset(ir,mod_group=="MRI"),
  US=subset(ir,mod_group=="CEUS/US"),Combined=subset(ir,mod_group=="CT/MRI"),
  CTMRI=subset(im,mod_pair=="CT vs MRI"),CTUS=subset(im,mod_pair=="CT vs CEUS/US"),
  CTCEUS=subset(im,mod_pair=="CT vs CEUS/US" & (std_modality_1=="CEUS"|std_modality_2=="CEUS")))
for (k in names(groups)) {
  d <- groups[[k]]; stopifnot(!anyDuplicated(d$study_id))
  m <- rma(yi,vi,data=d,method="REML")
  out[[length(out)+1L]] <- data.frame(analysis=paste("OPES",k),effects=m$k,clusters=m$k,studies=m$k,
    estimate=tanh(coef(m)[1]),ci_lower=tanh(m$ci.lb),ci_upper=tanh(m$ci.ub),pi_lower=NA_real_,pi_upper=NA_real_)
}
dest <- file.path(root,"revision/analysis/results")
write.csv(do.call(rbind,out),file.path(dest,"R2_full_table_additional_models.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(list(ir_res,im_res),function(x)
  data.frame(total=x$I2_total,between=x$I2_between,within=x$I2_within))) |>
  transform(stratum=c("inter-reader","inter-modality")),file.path(dest,"primary_I2_percent.csv"),row.names=FALSE)
