## Run from any directory. Read live R2 inputs; no submission/analysis writes.
arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
w <- dirname(normalizePath(arg)); root <- normalizePath(file.path(w,"../../.."))
.libPaths(c(file.path(root,"revision/R_lib"), .libPaths()))
suppressPackageStartupMessages({library(metafor); library(dplyr); library(readr)})
data_dir <- file.path(root,"meta_analysis"); out <- file.path(w,"data")
dir.create(out,showWarnings=FALSE)
## Evaluate ONLY the named input/preparation assignments from the approved
## R2 script. Do not execute its file writers or refit unrelated analyses.
needed <- c("comp","stud","qarel","pediatric_ids","dat_all","dat","ir","im","prepare_opes")
used <- character()
for(e in parse(file.path(root,"revision/analysis/revision_analyses.R"))) {
 if(is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) &&
    as.character(e[[2]]) %in% needed) { eval(e); used <- c(used,as.character(e[[2]])) }
}
stopifnot(identical(used,needed))
summary <- read_csv(file.path(root,"revision/analysis/results/revision_analysis_summary.csv"),show_col_types=FALSE)
rows <- boxes <- tests <- list()
for(s in c("inter-reader","inter-modality")) {
 d <- prepare_opes(if(s=="inter-reader") ir else im)
 stopifnot(nrow(d)==if(s=="inter-reader") 46 else 30,!anyDuplicated(d$study_id))
 if(s=="inter-reader") stopifnot(!("Graumann_2015" %in% d$study_id),
   identical(d$opes_kappa[d$study_id=="He_2026"],0.77))
 d$grp <- if(s=="inter-reader") {
   stopifnot(all(d$std_modality_1[d$mod_group=="other"]=="CT_MRI"))
   unname(c(CT="CT",MRI="MRI",CEUS_US="CEUS/US",other="Combined CT/MRI")[d$mod_group])
 } else unname(c(CT_MRI="CT vs MRI",CT_CEUS_US="CT vs US-based")[d$mod_pair])
 stopifnot(!anyNA(d$grp))
 model <- rma(opes_yi,opes_vi,data=d,method="REML")
 eg <- regtest(model,model="lm"); bg <- ranktest(model); tf <- trimfill(model)
 vals <- function(m) tanh(c(m$b[[1]], m$ci.lb, m$ci.ub))
 pooled <- vals(model); adjusted <- vals(tf)
 tests[[s]] <- tibble(stratum=s,k=nrow(d),egger_t=as.numeric(eg$zval),egger_df=eg$ddf,egger_p=eg$pval,
   begg_kendall_tau=bg$tau,begg_p=bg$pval,trimfill_imputed=tf$k0,trimfill_side=tf$side,
   opes_kappa=pooled[1],ci_lower=pooled[2],ci_upper=pooled[3],trimfill_kappa=adjusted[1],
   trimfill_ci_lower=adjusted[2],trimfill_ci_upper=adjusted[3])
 stopifnot(round(eg$pval,3)==if(s=="inter-reader") .212 else .309)
 if(s=="inter-reader") stopifnot(identical(round(pooled,3),c(.687,.613,.749)))
 one <- function(group,kind,k,v) tibble(stratum=s,group=group,type=kind,k=k,estimate=v[1],lower=v[2],upper=v[3])
 primary <- summary |> filter(analysis_family=="primary", startsWith(tolower(analysis),s))
 stopifnot(nrow(primary)==1)
 boxes[[s]] <- bind_rows(bind_rows(unname(lapply(split(d,d$grp),function(g)
    one(g$grp[1],"subtotal",nrow(g),vals(rma(opes_yi,opes_vi,data=g,method="REML")))))),
    one("Overall","OPES",nrow(d),pooled),
    one("Overall","RVE+CR2",primary$effects,c(primary$estimate,primary$ci_lower,primary$ci_upper)))
 rows[[s]] <- d |> transmute(stratum=s,study_id,comparison_id,group=grp,year,n=nc,
   yi=opes_yi,vi=opes_vi,estimate=tanh(opes_yi),lower=tanh(opes_yi-qnorm(.975)*sqrt(opes_vi)),
   upper=tanh(opes_yi+qnorm(.975)*sqrt(opes_vi)),weight=as.numeric(weights(model)),
   marker=case_when(s=="inter-modality" & (std_modality_1=="DECT"|std_modality_2=="DECT") ~ "a",
                    s=="inter-modality" & (std_modality_1=="CT_MRI"|std_modality_2=="CT_MRI") ~ "b", TRUE ~ ""))
}
write_csv(bind_rows(rows),file.path(out,"forest_rows.csv"))
write_csv(bind_rows(boxes),file.path(out,"forest_models.csv"))
write_csv(bind_rows(tests),file.path(out,"publication_bias_summary.csv"))
