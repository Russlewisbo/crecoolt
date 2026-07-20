## CRECOOLT — Novel vs. Older Regimens for CRE infection in liver transplant recipients
## Consolidated, reproducible analysis. Run:  source("analysis.R")
## Data:    CRECOOLT_onlyinfections.sav  (N = 405 CRE-infection cohort, 1 row/patient)
## Outputs: figure1_mortality_iptw.*, figure_3panel.*, figure_forest_90d.*,
##          loveplot_iptw.*, loveplot_era.*   (each as .png/.pdf/.svg)
## See STATUS.md for the narrative summary of findings.

suppressPackageStartupMessages({
  library(tidyverse); library(haven); library(survival); library(patchwork)
})

jama  <- c(Old = "#8C8C8C", Novel = "#00468B")   # JAMA grey / blue
save3 <- function(p, stem, w, h) for (e in c("png","pdf","svg"))
  ggsave(paste0(stem, ".", e), p, width = w, height = h, dpi = 300, bg = "white")

## ------------------------------------------------------------------ ##
## 1. Load and build the analysis frame                               ##
## ------------------------------------------------------------------ ##
new <- read_sav("CRECOOLT_onlyinfections.sav")

cr <- new |>
  transmute(
    record_id,
    regimen = factor(if_else(as.numeric(new_drug) == 1, "Novel", "Old"), levels = c("Old","Novel")),
    novel = as.integer(as.numeric(new_drug) == 1),
    study_period = as.numeric(study_period),
    age = as.numeric(age), meld = as.numeric(meld_score), sofa = as.numeric(sofa),
    kpc = as.numeric(mec_of_carbapenem_resi___1), bsi = as.numeric(infection_source___1),
    multisite = as.numeric(multisite_colonization),
    pgnf = as.numeric(post_olt_compli___4), crrt = as.numeric(post_olt_compli___2),
    combo = as.numeric(combination_treat),
    rec = as.numeric(cre_recurrence), res = as.numeric(resistance),
    t_rec = as.numeric(CREinf_to_relapse), t_res = as.numeric(CREinf_to_resistance),
    death = as.numeric(death), death90 = as.numeric(death_90d),
    t_death = as.numeric(CREinf_to_death),
    inf = as.Date(cre_infection_date), dc = as.Date(discharge_date),
    cl = as.Date(clearence_date), en = as.Date(end_date)
  ) |>
  mutate(
    fu = pmax(as.numeric(dc - inf), as.numeric(cl - inf), as.numeric(en - inf), na.rm = TRUE),
    last_valid = if_else(is.finite(fu) & fu >= 0, fu, NA_real_)   # 54 invalid discharge dates dropped
  )

## competing-risks (time, status): 1 = event of interest, 2 = death, 0 = censored
make_cr <- function(d, flag, tvar) d |>
  mutate(status = case_when(.data[[flag]] == 1 ~ 1L, death == 1 ~ 2L, TRUE ~ 0L),
         time = case_when(.data[[flag]] == 1 ~ .data[[tvar]], death == 1 ~ t_death, TRUE ~ last_valid),
         time = if_else(time == 0, 0.5, time)) |>
  filter(!is.na(time), is.finite(time), time > 0)

smd <- function(x, g, w = rep(1, length(x))) {
  m1 <- weighted.mean(x[g==1], w[g==1]); m0 <- weighted.mean(x[g==0], w[g==0])
  v1 <- weighted.mean((x[g==1]-m1)^2, w[g==1]); v0 <- weighted.mean((x[g==0]-m0)^2, w[g==0])
  (m1 - m0) / sqrt((v1 + v0) / 2)
}
fmt <- function(m, term = "regimenNovel") { s <- summary(m)
  c <- s$conf.int[term, c("exp(coef)","lower .95","upper .95")]
  sprintf("%.2f (%.2f-%.2f), p=%.3f", c[1], c[2], c[3], s$coefficients[term, ncol(s$coefficients)]) }

## ------------------------------------------------------------------ ##
## 2. Recurrence / resistance — competing risks (death competing)     ##
## ------------------------------------------------------------------ ##
d_rec <- make_cr(cr, "rec", "t_rec"); d_res <- make_cr(cr, "res", "t_res")

fg_full <- function(d, label) { d$ev <- factor(d$status, 0:2, c("censor", label, "death"))
  fg <- finegray(Surv(time, ev) ~ regimen, data = d[,c("time","ev","regimen")], etype = label)
  coxph(Surv(fgstart, fgstop, fgstatus) ~ regimen, weights = fgwt, data = fg) }

## era-restricted (2018-2025) IPTW Fine-Gray — the valid causal contrast
era2 <- cr |> filter(study_period == 2, !is.na(age))
era2$ps <- predict(glm(novel ~ age + meld + multisite + sofa, binomial, era2), type = "response")
pT <- mean(era2$novel); era2$iptw <- ifelse(era2$novel == 1, pT/era2$ps, (1-pT)/(1-era2$ps))
era_fg <- function(d, label) { d$ev <- factor(d$status, 0:2, c("censor", label, "death"))
  fg <- finegray(Surv(time, ev) ~ ., data = d[,c("time","ev","regimen","iptw","record_id")], etype = label)
  coxph(Surv(fgstart, fgstop, fgstatus) ~ regimen + cluster(record_id),
        weights = fg$fgwt * fg$iptw, data = fg, robust = TRUE) }
e_rec <- make_cr(era2, "rec", "t_rec"); e_res <- make_cr(era2, "res", "t_res")

## ------------------------------------------------------------------ ##
## 3. Mortality — 9-covariate IPTW (ATE), 90-day, from CRE infection  ##
## ------------------------------------------------------------------ ##
d2 <- cr |>
  mutate(time = case_when(death90 == 1 & !is.na(t_death) ~ pmin(t_death, 90),
                          TRUE ~ pmin(coalesce(last_valid, 90), 90)),
         time = if_else(time <= 0, 0.5, time),
         event = if_else(death90 == 1 & (is.na(t_death) | t_death <= 90), 1L, 0L)) |>
  drop_na(study_period, age, meld, sofa, kpc, bsi, multisite, pgnf, crrt)
d2$ps <- predict(glm(novel ~ study_period + age + meld + sofa + kpc + bsi + multisite + pgnf + crrt,
                     binomial, d2), type = "response")
d2$w    <- ifelse(d2$novel == 1, 1/d2$ps, 1/(1 - d2$ps))                 # ATE
d2$watt <- ifelse(d2$novel == 1, 1, d2$ps/(1 - d2$ps))                   # ATT (for diagnostics)
cox90 <- coxph(Surv(time, event) ~ regimen, d2, weights = w, robust = TRUE)

## ------------------------------------------------------------------ ##
## 4. Key results                                                     ##
## ------------------------------------------------------------------ ##
cat("\n================ KEY RESULTS ================\n")
cat("Cohort N =", nrow(new), "|", sum(cr$regimen=="Old"), "old /", sum(cr$regimen=="Novel"), "novel\n\n")
cat("MORTALITY (90-day, IPTW ATE):        HR", fmt(cox90), "\n")
cat("RECURRENCE full-cohort Fine-Gray:   sHR", fmt(fg_full(d_rec,"recurrence")), "\n")
cat("RESISTANCE full-cohort Fine-Gray:   sHR", fmt(fg_full(d_res,"resistance")), "\n")
cat("RECURRENCE 2018-25 IPTW Fine-Gray:  sHR", fmt(era_fg(e_rec,"recurrence")), "\n")
cat("RESISTANCE 2018-25 IPTW Fine-Gray:  sHR", fmt(era_fg(e_res,"resistance")), "\n")
cat("=============================================\n\n")

## ------------------------------------------------------------------ ##
## 5. Figure 1 — cumulative mortality (IPTW, 95% CI, number-at-risk)  ##
## ------------------------------------------------------------------ ##
sf9 <- survfit(Surv(time, event) ~ regimen, data = d2, weights = w, robust = TRUE)
s9  <- summary(sf9, times = 0:90, extend = TRUE)
curve_ci <- tibble(time = s9$time, ci = 1 - s9$surv, lower = 1 - s9$upper, upper = 1 - s9$lower,
                   regimen = factor(sub("regimen=", "", as.character(s9$strata)), levels = c("Old","Novel")))
lab  <- curve_ci |> group_by(regimen) |> filter(time == max(time)) |> ungroup()
tms  <- c(0,15,30,45,60,75,90); xexp <- expansion(add = c(7,4))
risk <- expand_grid(regimen = c("Old","Novel"), t = tms) |>
  mutate(n = map2_int(regimen, t, ~ sum(d2$regimen == .x & d2$time >= .y)),
         regimen = factor(regimen, levels = c("Old","Novel")))
mA <- ggplot(curve_ci, aes(time, ci, color = regimen, fill = regimen)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .18, colour = NA) + geom_step(linewidth = 1) +
  geom_text(data = lab, aes(label = regimen), hjust = -.05, size = 4, fontface = "bold") +
  annotate("text", x = 2, y = .74, hjust = 0, size = 3.8,
           label = paste0("IPTW-adjusted HR, ", sub(", p.*","",fmt(cox90)))) +
  scale_color_manual(values = jama) + scale_fill_manual(values = jama) +
  scale_y_continuous("Cumulative mortality, %", labels = \(x) x*100, limits = c(0,.80), expand = c(0,0)) +
  scale_x_continuous("Days since CRE infection", breaks = tms, limits = c(0,104), expand = xexp) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.line = element_line(linewidth = .4), plot.margin = margin(8,14,4,8))
mR <- ggplot(risk, aes(t, regimen, label = n, color = regimen)) + geom_text(size = 3.6) +
  scale_color_manual(values = jama) + scale_x_continuous(limits = c(0,104), breaks = tms, expand = xexp) +
  scale_y_discrete(limits = c("Novel","Old")) + labs(title = "No. at risk") + theme_void(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0),
        axis.text.y = element_text(color = jama[c("Novel","Old")], hjust = 1, size = 10), plot.margin = margin(0,14,4,8))
fig1 <- mA / mR + plot_layout(heights = c(6,1)); save3(fig1, "figure1_mortality_iptw", 7, 5.2)

## ------------------------------------------------------------------ ##
## 6. 3-panel figure (mortality | recurrence | resistance)            ##
## ------------------------------------------------------------------ ##
wcif <- function(d, label, nice) { d$ev <- factor(d$status, 0:2, c("censor", label, "death"))
  sf <- survfit(Surv(time, ev) ~ regimen, data = d, weights = iptw); k <- which(sf$states == label)
  s <- summary(sf, times = 0:90, extend = TRUE)
  tibble(time = s$time, est = s$pstate[,k], lower = s$lower[,k], upper = s$upper[,k],
         regimen = sub("regimen=", "", as.character(s$strata)), endpoint = nice) }
band <- bind_rows(wcif(e_rec,"recurrence","CRE recurrence"), wcif(e_res,"resistance","Resistance emergence")) |>
  mutate(regimen = factor(regimen, levels = c("Old","Novel")))
bt <- theme_classic(base_size = 12) + theme(legend.position = "top", axis.line = element_line(linewidth = .4),
        plot.title = element_text(face = "bold", size = 12), plot.margin = margin(6,12,4,6))
xs <- scale_x_continuous("Days since CRE infection", breaks = c(0,30,60,90), limits = c(0,90), expand = c(0,0))
pA <- mA + labs(title = "All-cause mortality") + bt + theme(legend.position = "none")
mk <- function(ep, ttl, m) ggplot(filter(band, endpoint == ep), aes(time, est, color = regimen, fill = regimen)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .18, colour = NA) + geom_step(linewidth = 1) +
  annotate("text", x = 2, y = .48, hjust = 0, size = 3.3, lineheight = 1,
           label = paste0("Fine-Gray sHR, ", sub(", p.*","",fmt(m)))) +
  scale_color_manual("Regimen", values = jama) + scale_fill_manual("Regimen", values = jama) +
  scale_y_continuous("Cumulative incidence, %", labels = \(x) x*100, limits = c(0,.52), expand = c(0,0)) +
  xs + labs(title = ttl) + bt
pB <- mk("CRE recurrence", "CRE recurrence", era_fg(e_rec,"recurrence"))
pC <- mk("Resistance emergence", "Resistance emergence", era_fg(e_res,"resistance"))
fig3 <- (pA | pB | pC) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "top"); save3(fig3, "figure_3panel", 12, 4.6)

## ------------------------------------------------------------------ ##
## 7. Love plots — full cohort (unbalanceable) & era (balanced)       ##
## ------------------------------------------------------------------ ##
love_plot <- function(dat, w, covs, labels, caption) {
  tibble(Covariate = labels,
         Before = map_dbl(names(covs), ~smd(dat[[.x]], dat$novel)),
         After  = map_dbl(names(covs), ~smd(dat[[.x]], dat$novel, w))) |>
    pivot_longer(c(Before, After), names_to = "phase", values_to = "smd") |>
    mutate(Covariate = fct_reorder(Covariate, abs(smd) * (phase == "Before")),
           phase = factor(phase, levels = c("Before","After"))) |>
    ggplot(aes(smd, Covariate, color = phase, shape = phase)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = .3) +
    geom_vline(xintercept = c(-.1,.1), linetype = "dashed", color = "grey60", linewidth = .3) +
    geom_point(size = 3) +
    scale_color_manual(NULL, values = c(Before = "#8C8C8C", After = "#00468B")) +
    scale_shape_manual(NULL, values = c(Before = 16, After = 17)) +
    labs(x = "Standardized mean difference", y = NULL, caption = caption) +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
}
covs9 <- c(study_period=1,multisite=1,meld=1,kpc=1,bsi=1,sofa=1,age=1,pgnf=1,crrt=1)
lab9  <- c("Study period","Multisite colonisation","MELD score","KPC carbapenemase",
           "Bloodstream infection","SOFA score","Age (years)","Primary graft non-function","Continuous RRT")
save3(love_plot(d2, d2$w, covs9, lab9, "Full cohort: study period not balanceable (positivity violation)"),
      "loveplot_iptw", 7.5, 4.6)
covs4 <- c(age=1,meld=1,multisite=1,sofa=1)
lab4  <- c("Age (years)","MELD score","Multisite colonisation","SOFA score")
save3(love_plot(era2, era2$iptw, covs4, lab4, "Restricted to 2018-2025 era: all covariates balanced"),
      "loveplot_era", 7, 4)

## ------------------------------------------------------------------ ##
## 8. Subgroup forest plot (90-day mortality, IPTW)                   ##
## ------------------------------------------------------------------ ##
fit_sub <- function(d) { s <- summary(coxph(Surv(time, event) ~ regimen, d, weights = w, robust = TRUE))
  c(s$conf.int["regimenNovel",1], s$conf.int["regimenNovel",3], s$conf.int["regimenNovel",4]) }
pint <- function(v) { m <- coxph(as.formula(paste0("Surv(time,event)~regimen*", v)), d2, weights = w, robust = TRUE)
  unname(tail(summary(m)$coefficients[,"Pr(>|z|)"], 1)) }
subs <- tribble(~group, ~label, ~flt,
  "Overall","Overall",quote(TRUE),
  "KPC status","KPC",quote(kpc==1),"KPC status","Non-KPC",quote(kpc==0),
  "Infection source","BSI",quote(bsi==1),"Infection source","Non-BSI",quote(bsi==0),
  "SOFA score","SOFA \u22647",quote(sofa<=7),"SOFA score","SOFA >7",quote(sofa>7),
  "Therapy type","Combination",quote(combo==1),"Therapy type","Monotherapy",quote(combo==0))
forest <- subs |> mutate(r = map(flt, \(e){ d <- d2[eval(e, d2), ]; c(n = nrow(d), fit_sub(d)) }),
  n = map_dbl(r,1), hr = map_dbl(r,2), lo = map_dbl(r,3), hi = map_dbl(r,4)) |> select(group,label,n,hr,lo,hi)
pv <- c("KPC status"=pint("kpc"),"Infection source"=pint("bsi"),
        "SOFA score"=pint("I(sofa>7)"),"Therapy type"=pint("combo"))
dr <- tribble(~row,~type,"Overall","est","","sp","KPC status","hdr","KPC","est","Non-KPC","est","","sp",
  "Infection source","hdr","BSI","est","Non-BSI","est","","sp","SOFA score","hdr","SOFA \u22647","est","SOFA >7","est",
  "","sp","Therapy type","hdr","Combination","est","Monotherapy","est") |>
  mutate(y = n() - row_number() + 1) |> left_join(select(forest, row = label, hr, lo, hi, n), by = "row") |>
  mutate(pint = setNames(unname(pv), c("KPC status","Infection source","SOFA score","Therapy type"))[row],
         lab = case_when(type=="hdr"~row, row=="Overall"~row, type=="est"~paste0("   ",row), TRUE~row),
         hrt = if_else(type=="est", sprintf("%.2f (%.2f-%.2f)",hr,lo,hi), NA),
         nt = if_else(type=="est", as.character(n), NA), pt = if_else(!is.na(pint), sprintf("%.2f",pint), NA))
tp <- max(dr$y) + 1.3
ptxt <- ggplot(dr, aes(y = y)) +
  geom_text(aes(x=0,label=lab,fontface=ifelse(type=="hdr","bold","plain")), hjust=0, size=3.3, na.rm=TRUE) +
  geom_text(aes(x=4.0,label=nt), hjust=1, size=3.1, na.rm=TRUE) +
  geom_text(aes(x=8.6,label=hrt), hjust=1, size=3.1, na.rm=TRUE) +
  geom_text(aes(x=9.9,label=pt), hjust=1, size=3.0, fontface="italic", color="grey35", na.rm=TRUE) +
  annotate("text",x=0,y=tp,label="Subgroup",hjust=0,fontface="bold",size=3.2) +
  annotate("text",x=4.0,y=tp,label="No.",hjust=1,fontface="bold",size=3.1) +
  annotate("text",x=8.6,y=tp,label="HR (95% CI)",hjust=1,fontface="bold",size=3.1) +
  annotate("text",x=9.9,y=tp,label="P-int",hjust=1,fontface="bold",size=3.0) +
  scale_x_continuous(limits=c(0,10.1)) + scale_y_continuous(limits=c(.5,tp+.5)) + theme_void()
pfor <- ggplot(filter(dr,type=="est"), aes(hr,y)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60", linewidth=.3) +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, color="#00468B", linewidth=.6, orientation="y") +
  geom_point(color="#00468B", size=2.3, shape=15) +
  annotate("text",x=.34,y=tp,label="Favors novel",size=2.7,color="grey40") +
  annotate("text",x=2.6,y=tp,label="Favors old",size=2.7,color="grey40") +
  scale_x_log10("Hazard ratio (95% CI)", breaks=c(.25,.5,1,2), limits=c(.2,3.6)) +
  scale_y_continuous(limits=c(.5,tp+.5)) + theme_classic(base_size=11) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.line.y=element_blank())
save3((ptxt | pfor) + plot_layout(widths = c(2.1,1)), "figure_forest_90d", 10, 4.6)

cat("All figures written (png/pdf/svg).\n")
