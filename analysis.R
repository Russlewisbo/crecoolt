## CRECOOLT — Novel vs. older regimens for CRE infection in liver transplant recipients
## Consolidated, reproducible analysis (2018-2024 cohort). Run:  source("analysis.R")
## Data:    CRECOOLT_onlyinfections.sav (N=405); analyses restricted to the 2018-2025 era (N=182).
## Design:  PRIMARY  = 30-day all-cause mortality from CRE infection (IPTW-weighted Cox)
##          SECONDARY = 90-day CRE recurrence & resistance emergence (IPTW Fine-Gray, death competing)
##          PS covariates: age, SOFA, BSI, multisite colonization, MELD (stabilized ATE weights)
## Outputs: figure1_mortality_iptw.*, figure2_recurrence_resistance.*, figure_forest_30d.*,
##          loveplot_era.*  (each png/pdf/svg). See STATUS.md for narrative.

suppressPackageStartupMessages({ library(tidyverse); library(haven); library(survival); library(patchwork) })
jama  <- c(Old = "#8C8C8C", Novel = "#00468B")
save3 <- function(p, stem, w, h) for (e in c("png","pdf","svg"))
  ggsave(paste0(stem, ".", e), p, width = w, height = h, dpi = 300, bg = "white")
fmt <- function(m, term = "regimenNovel") { s <- summary(m)
  c <- s$conf.int[term, c("exp(coef)","lower .95","upper .95")]
  sprintf("%.2f (%.2f-%.2f), p=%.3f", c[1], c[2], c[3], s$coefficients[term, ncol(s$coefficients)]) }

## ---- 1. Build the 2018-2025 era analysis frame ----
new <- read_sav("CRECOOLT_onlyinfections.sav")
dat <- new |>
  filter(as.numeric(study_period) == 2) |>
  transmute(
    record_id,
    regimen = factor(if_else(as.numeric(new_drug) == 1, "Novel", "Old"), levels = c("Old","Novel")),
    novel = as.integer(as.numeric(new_drug) == 1),
    age = as.numeric(age), meld = as.numeric(meld_score), sofa = as.numeric(sofa),
    multisite = as.numeric(multisite_colonization), bsi = as.numeric(infection_source___1),
    kpc = as.numeric(mec_of_carbapenem_resi___1), combo = as.numeric(combination_treat),
    rec = as.numeric(cre_recurrence), res = as.numeric(resistance),
    t_rec = as.numeric(CREinf_to_relapse), t_res = as.numeric(CREinf_to_resistance),
    death = as.numeric(death), death30 = as.numeric(death_30d),
    inf = as.Date(cre_infection_date), dth = as.Date(death_date),
    dc = as.Date(discharge_date), cl = as.Date(clearence_date), en = as.Date(end_date)
  ) |>
  mutate(
    t_death = as.numeric(dth - inf),                       # date-derived competing/mortality time
    fu = pmax(as.numeric(dc-inf), as.numeric(cl-inf), as.numeric(en-inf), na.rm = TRUE),
    last_valid = if_else(is.finite(fu) & fu >= 0, fu, NA_real_)
  ) |>
  drop_na(age, sofa, bsi, multisite, meld)                 # complete cases on PS covariates

# time from infection to first regimen start (for immortal-time / delay sensitivity)
dat <- dat |>
  left_join(new |> transmute(record_id, t_treat = as.numeric(as.Date(start_date) - as.Date(cre_infection_date))),
            by = "record_id") |>
  mutate(t_treat0 = if_else(is.na(t_treat) | t_treat < 0, 0, t_treat))   # empiric/negative delays -> 0

## ---- 2. Propensity score & stabilized IPTW (ATE) ----
dat$ps <- predict(glm(novel ~ age + sofa + bsi + multisite + meld, binomial, dat), type = "response")
pT <- mean(dat$novel)
dat$iptw <- ifelse(dat$novel == 1, pT / dat$ps, (1 - pT) / (1 - dat$ps))

## ---- 3. Time/status builders ----
# 30-day mortality (primary)
dat <- dat |> mutate(
  mtime = case_when(death30 == 1 & !is.na(t_death) ~ pmin(t_death, 30), TRUE ~ pmin(coalesce(last_valid, 30), 30)),
  mtime = if_else(mtime <= 0, 0.5, mtime),
  mevent = if_else(death30 == 1 & (is.na(t_death) | t_death <= 30), 1L, 0L))
# competing-risks (time,status): 1=event, 2=death, 0=censored (secondary, 90-day)
make_cr <- function(d, flag, tvar) d |>
  mutate(status = case_when(.data[[flag]] == 1 ~ 1L, death == 1 ~ 2L, TRUE ~ 0L),
         time = case_when(.data[[flag]] == 1 ~ .data[[tvar]], death == 1 ~ t_death, TRUE ~ last_valid),
         time = if_else(time == 0, 0.5, time)) |>
  filter(!is.na(time), is.finite(time), time > 0)
d_rec <- make_cr(dat, "rec", "t_rec"); d_res <- make_cr(dat, "res", "t_res")

## ---- 4. Models ----
cox30 <- coxph(Surv(mtime, mevent) ~ regimen, dat, weights = iptw, robust = TRUE)   # primary
fg <- function(d, label) { d$ev <- factor(d$status, 0:2, c("censor", label, "death"))
  f <- finegray(Surv(time, ev) ~ ., data = d[,c("time","ev","regimen","iptw","record_id")], etype = label)
  coxph(Surv(fgstart, fgstop, fgstatus) ~ regimen + cluster(record_id), weights = f$fgwt*f$iptw, data = f, robust = TRUE) }
sh_rec <- fg(d_rec, "recurrence"); sh_res <- fg(d_res, "resistance")

cat("\n=============== KEY RESULTS (2018-2025 era, N =", nrow(dat), ") ===============\n")
cat("30-day deaths:", sum(dat$mevent), "| crude old",
    scales::percent(mean(dat$death30[dat$regimen=="Old"]),.1), "vs novel",
    scales::percent(mean(dat$death30[dat$regimen=="Novel"]),.1), "\n")
cat("PRIMARY  30-day mortality, IPTW Cox HR:  ", fmt(cox30), "\n")
cat("SECONDARY recurrence,  IPTW Fine-Gray sHR:", fmt(sh_rec), "\n")
cat("SECONDARY resistance,  IPTW Fine-Gray sHR:", fmt(sh_res), "\n")
cat("========================================================================\n\n")

## ---- 5. Figure 1: 30-day cumulative mortality ----
sf <- survfit(Surv(mtime, mevent) ~ regimen, data = dat, weights = iptw, robust = TRUE)
s <- summary(sf, times = 0:30, extend = TRUE)
cm <- tibble(time = s$time, ci = 1 - s$surv, lower = 1 - s$upper, upper = 1 - s$lower,
             regimen = factor(sub("regimen=", "", as.character(s$strata)), levels = c("Old","Novel")))
lab <- cm |> group_by(regimen) |> filter(time == max(time)) |> ungroup()
tms <- c(0,5,10,15,20,25,30); xexp <- expansion(add = c(2.5,1.5))
risk <- expand_grid(regimen = c("Old","Novel"), t = tms) |>
  mutate(n = map2_int(regimen, t, ~ sum(dat$regimen == .x & dat$mtime >= .y)), regimen = factor(regimen, levels = c("Old","Novel")))
mA <- ggplot(cm, aes(time, ci, color = regimen, fill = regimen)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .18, colour = NA) + geom_step(linewidth = 1) +
  geom_text(data = lab, aes(label = regimen), hjust = -.05, size = 4, fontface = "bold") +
  annotate("text", x = .5, y = .52, hjust = 0, size = 3.8, label = paste0("IPTW-adjusted HR, ", sub(", p.*","",fmt(cox30)))) +
  scale_color_manual(values = jama) + scale_fill_manual(values = jama) +
  scale_y_continuous("Cumulative mortality, %", labels = \(x) x*100, limits = c(0,.56), expand = c(0,0)) +
  scale_x_continuous("Days since CRE infection", breaks = tms, limits = c(0,35), expand = xexp) +
  labs(caption = "30-day all-cause mortality; IPTW-weighted, 2018\u20132025 era (N=182)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.line = element_line(linewidth = .4),
        plot.caption = element_text(hjust = 0, color = "grey40"), plot.margin = margin(8,14,2,8))
mR <- ggplot(risk, aes(t, regimen, label = n, color = regimen)) + geom_text(size = 3.6) +
  scale_color_manual(values = jama) + scale_x_continuous(limits = c(0,35), breaks = tms, expand = xexp) +
  scale_y_discrete(limits = c("Novel","Old")) + labs(title = "No. at risk") + theme_void(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0),
        axis.text.y = element_text(color = jama[c("Novel","Old")], hjust = 1, size = 10), plot.margin = margin(0,14,4,8))
save3(mA / mR + plot_layout(heights = c(6,1)), "figure1_mortality_iptw", 7, 5.2)

## ---- 6. Figure 2: recurrence & resistance (90-day, competing risks) ----
wcif <- function(d, label, nice) { d$ev <- factor(d$status, 0:2, c("censor", label, "death"))
  sf <- survfit(Surv(time, ev) ~ regimen, data = d, weights = iptw); k <- which(sf$states == label)
  s <- summary(sf, times = 0:90, extend = TRUE)
  tibble(time = s$time, est = s$pstate[,k], lower = s$lower[,k], upper = s$upper[,k],
         regimen = sub("regimen=", "", as.character(s$strata)), endpoint = nice) }
band <- bind_rows(wcif(d_rec,"recurrence","CRE recurrence"), wcif(d_res,"resistance","Resistance emergence")) |>
  mutate(regimen = factor(regimen, levels = c("Old","Novel")))
ann <- tibble(endpoint = c("CRE recurrence","Resistance emergence"),
              lab = c(paste0("Fine-Gray sHR\n", sub(", p.*","",fmt(sh_rec))),
                      paste0("Fine-Gray sHR\n", sub(", p.*","",fmt(sh_res)))))
p2 <- ggplot(band, aes(time, est, color = regimen, fill = regimen)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .18, colour = NA) + geom_step(linewidth = 1) +
  geom_text(data = ann, aes(x = 3, y = .48, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.4, lineheight = .95) +
  facet_wrap(~ endpoint) + scale_color_manual("Regimen", values = jama) + scale_fill_manual("Regimen", values = jama) +
  scale_y_continuous("Cumulative incidence, %", labels = \(x) x*100, limits = c(0,.52), expand = c(0,0)) +
  scale_x_continuous("Days since CRE infection", breaks = c(0,30,60,90), limits = c(0,90), expand = expansion(mult = c(.03,.06))) +
  labs(caption = "IPTW-weighted, 2018\u20132025 era (N=182); death as competing risk") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", axis.line = element_line(linewidth = .4), panel.spacing.x = unit(2, "lines"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, color = "grey40", size = 9))
save3(p2, "figure2_recurrence_resistance", 8.7, 4.4)

## ---- 7. Love plot: 5-covariate PS balance ----
smd <- function(x, g, w = rep(1, length(x))) { m1<-weighted.mean(x[g==1],w[g==1]); m0<-weighted.mean(x[g==0],w[g==0])
  v1<-weighted.mean((x[g==1]-m1)^2,w[g==1]); v0<-weighted.mean((x[g==0]-m0)^2,w[g==0]); (m1-m0)/sqrt((v1+v0)/2) }
cv <- c(age="Age", sofa="SOFA score", bsi="Bloodstream infection", multisite="Multisite colonisation", meld="MELD score")
lp <- tibble(Covariate = cv, Before = map_dbl(names(cv), ~smd(dat[[.x]], dat$novel)),
             After = map_dbl(names(cv), ~smd(dat[[.x]], dat$novel, dat$iptw))) |>
  pivot_longer(c(Before,After), names_to="phase", values_to="smd") |>
  mutate(Covariate = fct_reorder(Covariate, abs(smd)*(phase=="Before")), phase = factor(phase, levels=c("Before","After"))) |>
  ggplot(aes(smd, Covariate, color = phase, shape = phase)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = .3) +
  geom_vline(xintercept = c(-.1,.1), linetype = "dashed", color = "grey60", linewidth = .3) + geom_point(size = 3) +
  scale_color_manual(NULL, values = c(Before="#8C8C8C", After="#00468B")) +
  scale_shape_manual(NULL, values = c(Before=16, After=17)) +
  labs(x = "Standardized mean difference", y = NULL, caption = "2018-2025 era (N=182); PS: age, SOFA, BSI, multisite colonisation, MELD") +
  theme_minimal(base_size = 12) + theme(legend.position = "top", plot.caption = element_text(hjust = 0, color = "grey40"))
save3(lp, "loveplot_era", 7, 3.6)

## ---- 8. Forest: 30-day mortality subgroups ----
fit_s <- function(d) { m <- tryCatch(coxph(Surv(mtime,mevent)~regimen, d, weights=iptw, robust=TRUE), error=\(e) NULL)
  if (is.null(m)) return(c(NA,NA,NA)); summary(m)$conf.int["regimenNovel", c("exp(coef)","lower .95","upper .95")] }
pint <- function(v) { m <- tryCatch(coxph(as.formula(paste0("Surv(mtime,mevent)~regimen*",v)), dat, weights=iptw, robust=TRUE), error=\(e) NULL)
  if (is.null(m)) return(NA_real_); unname(tail(summary(m)$coefficients[,"Pr(>|z|)"],1)) }
subs <- tribble(~group,~label,~flt,"Overall","Overall",quote(rep(TRUE,nrow(dat))),
  "KPC status","KPC",quote(kpc==1),"KPC status","Non-KPC",quote(kpc==0),
  "Infection source","BSI",quote(bsi==1),"Infection source","Non-BSI",quote(bsi==0),
  "SOFA score","SOFA \u22647",quote(sofa<=7),"SOFA score","SOFA >7",quote(sofa>7),
  "Therapy type","Combination",quote(combo==1),"Therapy type","Monotherapy",quote(combo==0))
fe <- subs |> mutate(r = map(flt, \(e){ d <- dat[eval(e, dat),]; c(n=nrow(d), fit_s(d)) }),
  n = map_dbl(r,1), hr = map_dbl(r,2), lo = map_dbl(r,3), hi = map_dbl(r,4)) |> select(group,label,n,hr,lo,hi)
pv <- c("KPC status"=pint("kpc"),"Infection source"=pint("bsi"),"SOFA score"=pint("I(sofa>7)"),"Therapy type"=pint("combo"))
dr <- tribble(~row,~type,"Overall","est","","sp","KPC status","hdr","KPC","est","Non-KPC","est","","sp",
  "Infection source","hdr","BSI","est","Non-BSI","est","","sp","SOFA score","hdr","SOFA \u22647","est","SOFA >7","est",
  "","sp","Therapy type","hdr","Combination","est","Monotherapy","est") |>
  mutate(y = n() - row_number() + 1) |> left_join(select(fe, row=label, hr, lo, hi, n), by="row") |>
  mutate(pint = setNames(unname(pv), c("KPC status","Infection source","SOFA score","Therapy type"))[row],
         lab = case_when(type=="hdr"~row, row=="Overall"~row, type=="est"~paste0("   ",row), TRUE~row),
         hrt = if_else(type=="est" & !is.na(hr), sprintf("%.2f (%.2f-%.2f)",hr,lo,hi), NA),
         nt = if_else(type=="est", as.character(n), NA), pt = if_else(!is.na(pint), sprintf("%.2f",pint), NA))
tp <- max(dr$y) + 1.3
pt_p <- ggplot(dr, aes(y=y)) +
  geom_text(aes(x=0,label=lab,fontface=ifelse(type=="hdr","bold","plain")),hjust=0,size=3.3,na.rm=TRUE) +
  geom_text(aes(x=4.0,label=nt),hjust=1,size=3.1,na.rm=TRUE) + geom_text(aes(x=8.6,label=hrt),hjust=1,size=3.1,na.rm=TRUE) +
  geom_text(aes(x=9.9,label=pt),hjust=1,size=3.0,fontface="italic",color="grey35",na.rm=TRUE) +
  annotate("text",x=0,y=tp,label="Subgroup",hjust=0,fontface="bold",size=3.2) +
  annotate("text",x=4.0,y=tp,label="No.",hjust=1,fontface="bold",size=3.1) +
  annotate("text",x=8.6,y=tp,label="HR (95% CI)",hjust=1,fontface="bold",size=3.1) +
  annotate("text",x=9.9,y=tp,label="P-int",hjust=1,fontface="bold",size=3.0) +
  scale_x_continuous(limits=c(0,10.1)) + scale_y_continuous(limits=c(.5,tp+.5)) + theme_void()
pf <- ggplot(filter(dr,type=="est"), aes(hr,y)) +
  geom_vline(xintercept=1,linetype="dashed",color="grey60",linewidth=.3) +
  geom_errorbar(aes(xmin=lo,xmax=hi),width=0,color="#00468B",linewidth=.6,orientation="y") +
  geom_point(color="#00468B",size=2.3,shape=15) +
  annotate("text",x=.28,y=tp,label="Favors novel",size=2.7,color="grey40") +
  annotate("text",x=3,y=tp,label="Favors old",size=2.7,color="grey40") +
  scale_x_log10("Hazard ratio (95% CI)",breaks=c(.25,.5,1,2),limits=c(.12,5)) +
  scale_y_continuous(limits=c(.5,tp+.5)) + theme_classic(base_size=11) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.line.y=element_blank())
save3((pt_p | pf) + plot_layout(widths=c(2.1,1)) +
        plot_annotation(caption="30-day mortality, IPTW-weighted, 2018-2025 era (N=182)",
                        theme=theme(plot.caption=element_text(hjust=0,color="grey40",size=9))),
      "figure_forest_30d", 10, 4.6)
## ---- 9. Time-to-treatment / immortal-time sensitivity ----
gethr <- function(m) summary(m)$conf.int["regimenNovel", c("exp(coef)","lower .95","upper .95")]
lm_keep <- function(d) filter(d, death == 0 | is.na(t_treat) | t_death >= t_treat)   # alive at treatment start
fg_s <- function(d, label, rhs = "regimen") {
  d$ev <- factor(d$status, 0:2, c("censor", label, "death"))
  vars <- unique(c("time","ev","iptw","record_id", all.vars(as.formula(paste("~", rhs)))))
  f <- finegray(Surv(time, ev) ~ ., data = d[, vars], etype = label)
  coxph(as.formula(paste0("Surv(fgstart,fgstop,fgstatus)~", rhs, "+cluster(record_id)")),
        weights = f$fgwt * f$iptw, data = f, robust = TRUE) }

cat("\n--- Immortal-time check: deaths before treatment start:",
    sum(dat$mevent == 1 & !is.na(dat$t_treat) & dat$t_death < dat$t_treat, na.rm = TRUE), "of", sum(dat$mevent), "---\n")

sens <- tibble::tribble(
  ~endpoint, ~analysis, ~m,
  "30-day mortality","Primary (IPTW)", cox30,
  "30-day mortality","+ time-to-treatment", coxph(Surv(mtime,mevent)~regimen+t_treat0, dat, weights=iptw, robust=TRUE),
  "30-day mortality","Landmark (alive at Tx)", coxph(Surv(mtime,mevent)~regimen, filter(dat, mevent==0 | is.na(t_treat) | t_death>=t_treat), weights=iptw, robust=TRUE),
  "Recurrence (90d)","Primary (IPTW)", sh_rec,
  "Recurrence (90d)","+ time-to-treatment", fg_s(d_rec,"recurrence","regimen+t_treat0"),
  "Recurrence (90d)","Landmark (alive at Tx)", fg_s(lm_keep(d_rec),"recurrence"),
  "Resistance (90d)","Primary (IPTW)", sh_res,
  "Resistance (90d)","+ time-to-treatment", fg_s(d_res,"resistance","regimen+t_treat0"),
  "Resistance (90d)","Landmark (alive at Tx)", fg_s(lm_keep(d_res),"resistance")) |>
  mutate(e = map(m, gethr), hr = map_dbl(e,1), lo = map_dbl(e,2), hi = map_dbl(e,3)) |>
  select(endpoint, analysis, hr, lo, hi)
cat("\n--- Time-to-treatment sensitivity ---\n")
print(sens |> mutate(txt = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)) |> select(endpoint, analysis, txt), n = Inf)

## ---- 10. Sensitivity figures ----
ttfig <- ggplot(dat, aes(regimen, t_treat, color = regimen, fill = regimen)) +
  geom_boxplot(width = .5, alpha = .2, outlier.shape = NA) +
  geom_jitter(width = .12, height = .2, size = 1, alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = .3) +
  scale_color_manual(values = jama) + scale_fill_manual(values = jama) +
  coord_cartesian(ylim = c(-5, 25)) +
  labs(x = NULL, y = "Days from CRE infection to regimen start",
       caption = "Dashed line = infection (culture) date; negative = empiric therapy started earlier") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.line = element_line(linewidth = .4),
        plot.caption = element_text(hjust = 0, color = "grey40", size = 9))
save3(ttfig, "figure_time_to_treatment", 5.5, 4)

sensf <- sens |>
  mutate(analysis = factor(analysis, levels = c("Landmark (alive at Tx)","+ time-to-treatment","Primary (IPTW)")),
         endpoint = factor(endpoint, levels = c("30-day mortality","Recurrence (90d)","Resistance (90d)")))
sfig <- ggplot(sensf, aes(hr, analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .18, color = "#00468B", linewidth = .6) +
  geom_point(color = "#00468B", size = 2.4, shape = 15) +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)), vjust = -.9, size = 2.9, color = "grey20") +
  facet_wrap(~ endpoint, ncol = 1, scales = "free_y") +
  scale_x_log10("HR (mortality) / subdistribution HR (recurrence, resistance), 95% CI", breaks = c(.25,.5,1,2,4), limits = c(.2,9)) +
  labs(y = NULL) + theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", color = NA), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
save3(sfig, "figure_sensitivity", 8, 5)

cat("All figures written (png/pdf/svg).\n")
