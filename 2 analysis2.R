################################################################################
##
##  Annual Ozone, Temperature, Mortality, Population data in Japan 
##  Chris Fook Sheng Ng, chrisng@m.u-tokyo.ac.jp
##  Matthew Lam
##
################################################################################

load(file="mort_tky_year.RData")
#mort_env <- mort_env[mort_env$year>1978,]

##### Summary Statistics #####

# visual
library(dplyr)
library(ggplot2)
library(patchwork)

# --- 1. Aggregate yearly ---
annual_summary <- mort_env %>%
  group_by(year) %>%
  summarise(
    Deaths      = sum(Deaths, na.rm = TRUE),
    PopulationM = sum(Population, na.rm = TRUE) / 1e6,                 # millions
    Ozone_ppb   = mean(mean_avg_daily_peak_ppm, na.rm = TRUE) * 1000,  # ppm → ppb
    Temp_C      = mean(daily_mean, na.rm = TRUE),                       # °C
    .groups = "drop"
  )

# --- 2. Individual plots with short titles + y-axis labels ---
p_deaths <- ggplot(annual_summary, aes(x = year, y = Deaths)) +
  geom_line(color = "red", linewidth = 1) +
  labs(title = "Mortality", y = "count", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

p_pop <- ggplot(annual_summary, aes(x = year, y = PopulationM)) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Population", y = "millions", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

p_o3 <- ggplot(annual_summary, aes(x = year, y = Ozone_ppb)) +
  geom_line(color = "purple", linewidth = 1) +
  labs(title = "Oxidant", y = "ppb", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

p_temp <- ggplot(annual_summary, aes(x = year, y = Temp_C)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Temperature", y = "°C", x = "Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

# Combine with main title
patch_plot <- (p_deaths / p_pop / p_o3 / p_temp) +
  plot_annotation(
    title = "Annual Mortality, Population, Oxidant, and Temperature in Tokyo (1976–2022)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  axis.title.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12)
                  )
  )
patch_plot

# Save to PDF
ggsave("tokyo_mort_ox.pdf", plot = patch_plot,
       width = 7.5, height = 10, units = "in")


##### Data Analysis #####
# --- 0) pick ozone metric & prep factors/offset/scaling ---
# Choose one ozone metric to avoid collinearity between ozone columns
ozone_var <- "mean_avg_daily_peak_ppm_lag1"  # the exposure one year before death

mort_env_mod <- mort_env |>
  mutate(
    # ensure factors are set (your data already has them, but make explicit)
    age_group = factor(age_group,
                       levels = c("0","1-4","5-14","15-24","25-44",
                                  "45-64","65-74","75-84","85-94","95+"),
                       ordered = FALSE),
    
    age_group5 = case_when(
      age_group %in% c("0", "1-4", "5-14", "15-24", "25-44", "45-64") ~ "0-64",
      age_group %in% c("65-74") ~ "65-74",
      age_group %in% c("75-84") ~ "75-84",
      age_group %in% c("85-94") ~ "85-94",
      age_group %in% c("95+") ~ "95+"
      ),
    age_group5 = factor(
      age_group5,
      levels = c("0-64","65-74","75-84","85-94","95+")),

    age_group7 = case_when(
      age_group %in% c("0", "1-4", "5-14") ~ "0-14",
      age_group %in% c("15-24", "25-44") ~ "15-44",
      age_group %in% c("45-64") ~ "45-64",
      age_group %in% c("65-74") ~ "65-74",
      age_group %in% c("75-84") ~ "75-84",
      age_group %in% c("85-94") ~ "85-94",
      age_group %in% c("95+") ~ "95+"
    ),
    age_group7 = factor(
      age_group7,
      levels = c("0-14","15-44","45-64","65-74","75-84","85-94","95+")),
    
        sex = factor(sex, levels = c("Female","Male")),
    # scale ozone to 0.01 ppm units (10 ppb) and center for stability
    ozone_10ppb = .data[[ozone_var]] * 100,
    ozone_10ppb_c = scale(ozone_10ppb, center = TRUE, scale = FALSE)[, 1],
    # center temperature too (use daily_mean to start)
    temp_c = scale(daily_max_lag1, center = TRUE, scale = FALSE)[, 1],
    # optional: keep to overlap years only, if you haven't already
    # year between 1976 and 2022 was your overlap assumption
    # year filter not applied here if mort_env already restricted
    offset_log_pop = log(Population)
  ) |>
  filter(!is.na(Deaths), !is.na(Population), Population > 0)

################################################################################
### Overall ozone effect
m_qp <- glm(
  Deaths ~ ozone_10ppb_c + sex + age_group7 + temp_c + poly(year,3) +
    offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)
summary(m_qp)

# Note: 
# Ozone coefficient is interpreted as the rate ratio per 10ppb increase 
# because we scaled by *100. (1ppm = 1000ppb)
# poly(year, 3) gives a smooth long-term trend in the GLM

# Extract coefficients and standard errors
coefs <- coef(m_qp)
ses   <- sqrt(diag(vcov(m_qp)))

# Compute 95% CI on the log scale
lower <- coefs - 1.96 * ses
upper <- coefs + 1.96 * ses

# Transform to RR scale
results <- data.frame(
  Term    = names(coefs),
  logRR   = coefs,
  SE      = ses,
  RR      = exp(coefs),
  LCL     = exp(lower),
  UCL     = exp(upper)
)

print(results) # check the results for ozone per 10ppb, RR and 95%CI
# RR = 1.0055 (95% CI: 0.9510, 1.0630)
# The RR is above 1, suggesting a positive association. However,
# the 95% CI contains 1 and the p-value = 0.847 > 0.05,
# suggesting the association is not statistically significant.
# Annual daytime ozone in the previous year (lag1) does not increase the risk 
# of annual all-cause mortality

# For high temperature (heat), represented by annual average daily max
# temperature, the RR (for temp_c) is 1.0124 (95% CI: 0.9781, 1.0479)
# The RR is > 1 showing a positive association of annual heat in the previous
# year (lag 1) with annual all-cause mortality rate. However, the 95%CI
# contains 1 and p = 0.484 > 0.05, suggesting the association is not
# statistical significance 


################################################################################
###  By sex
m_qp_sex <- glm(
  Deaths ~ ozone_10ppb_c*sex + age_group7 + temp_c + poly(year,3) +
    offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)

summary(m_qp_sex) # note the p-value for ozone_10ppb_c:sexMale is 0.0262 <0.05,
# which is statistically significant, suggesting a significant difference between sexes


co <- coef(m_qp_sex)
vc <- vcov(m_qp_sex)

# Female
b_f  <- co["ozone_10ppb_c"]
se_f <- sqrt(vc["ozone_10ppb_c", "ozone_10ppb_c"])
z_f  <- b_f / se_f
p_f  <- 2 * pnorm(-abs(z_f))
RR_f <- exp(b_f)
CI_f <- exp(b_f + c(-1.96, 1.96) * se_f)

# Male
b_m <- b_f + co["ozone_10ppb_c:sexMale"]
var_m <- vc["ozone_10ppb_c", "ozone_10ppb_c"] +
  vc["ozone_10ppb_c:sexMale", "ozone_10ppb_c:sexMale"] +
  2 * vc["ozone_10ppb_c", "ozone_10ppb_c:sexMale"]
se_m <- sqrt(var_m)
z_m  <- b_m / se_m
p_m  <- 2 * pnorm(-abs(z_m))
RR_m <- exp(b_m)
CI_m <- exp(b_m + c(-1.96, 1.96) * se_m)

# tidy table
data.frame(
  sex  = c("Female", "Male"),
  logRR = c(b_f, b_m),
  SE    = c(se_f, se_m),
  z     = c(z_f, z_m),
  pval  = c(p_f, p_m),
  RR    = c(RR_f, RR_m),
  LCL   = c(CI_f[1], CI_m[1]),
  UCL   = c(CI_f[2], CI_m[2])
)
# check the sex-specific results for ozone
# males show a positive association with RR >1, while females show a negative
# association with RR <1. However, the effect is only significant for male, 
# where its 95% CI does not contain 1, and p-value < 0.05.
# they are significantly different between the sexes.
# As noted above, the interaction term (ozone_10ppb_c:sexMale) shows a p-value < 0.05
# suggesting ozone effects differ by sex. 
# Results show male RR is 1.0323 (95%CI: 0.9720, 1.0963) which contains 1 and is
# not statistically significant
# Female RR is lower at 0.9717 (95%CI: 0.9124, 1.0348) which contains 1 and is 
# also not statistically significant

################################################################################
### By age (7 age groups)
m_qp_age7 <- glm(
  Deaths ~ ozone_10ppb_c * age_group7 + sex + temp_c + poly(year,3) +
    offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)

summary(m_qp_age7)

# --- manual RR per +10 ppb by 10 age groups (ref = "0") ---
coefs <- coef(m_qp_age7)     # <- your quasi-Poisson model object
vc    <- vcov(m_qp_age7)

base_var <- "ozone_10ppb_c"   # slope name for reference age group
age_fac  <- "age_group7"     # factor name used in the model
levels7 <- c("0-14","15-44","45-64","65-74","75-84","85-94","95+")

# helpers (safe getters)
getc <- function(x, n) if (n %in% names(x)) x[[n]] else NA_real_
getv <- function(M, r, c) if (r %in% rownames(M) && c %in% colnames(M)) M[r, c] else NA_real_

# containers
k   <- length(levels7)
est <- se <- numeric(k)

# reference group ("0")
est[1] <- getc(coefs, base_var)
se[1]  <- sqrt(getv(vc, base_var, base_var))

# other groups: add interaction terms base + int
for (i in 2:k) {
  lev <- levels7[i]
  int_name <- paste0(base_var, ":", age_fac, lev)  # e.g. "ozone_10ppb_c:age_group1065-74"
  
  beta_base <- getc(coefs, base_var)
  beta_int  <- getc(coefs, int_name)
  
  est[i] <- beta_base + beta_int
  
  var_i <- getv(vc, base_var, base_var) +
    getv(vc, int_name, int_name) +
    2 * getv(vc, base_var, int_name)
  
  se[i] <- sqrt(var_i)
}

# convert to RR and 95% CI
RR  <- exp(est)
LCL <- exp(est - 1.96 * se)
UCL <- exp(est + 1.96 * se)

# Wald test z and p-values
zval <- est / se
pval <- 2 * (1 - pnorm(abs(zval)))

# tidy table
age_rr7 <- data.frame(
  age_group = levels7,
  logRR     = est,
  SE        = se,
  RR        = RR,
  LCL       = LCL,
  UCL       = UCL,
  z         = zval,
  p         = pval
)

age_rr7
# check the age-specific results for ozone
# Show RR and 95% CI in a plot
library(dplyr)
library(ggplot2)

# keep the current data-frame order
age_order <- age_rr7$age_group

age_rr7 <- age_rr7 %>%
  mutate(age_group = factor(age_group, levels = age_order))

age_plot <- ggplot(age_rr7, aes(x = age_group, y = RR)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Age group", y = "Relative Risk (per 10 ppb ozone)") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(drop = FALSE)  # keep all levels even if missing

age_plot

# Save to PDF
ggsave("age_mort_ox.pdf", plot = age_plot,
       width = 6, height = 4, units = "in")


# Linear trend test
mort_env_mod <- mort_env_mod %>%
  mutate(age_ord = as.numeric(factor(age_group, 
                                     levels = c("0-14","15-44","45-64","65-74",
                                                "75-84","85-94","95+"))))

m_qp_trend <- glm(
  Deaths ~ ozone_10ppb_c * age_ord + sex + temp_c + poly(year,3) + 
    offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)

summary(m_qp_trend)

# However, linear trend test shows p-value 0.0997 >0.05, providing no evidence
# to support a linear trend of ozone risk across time.


################################################################################
### By age and sex 

m_qp_age7_sex <- glm(
  Deaths ~ ozone_10ppb_c * age_group7 * sex + temp_c + poly(year, 3) +
    offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)

summary(m_qp_age7_sex)

mort_env_mod <- mort_env_mod |>
  mutate(
    age_group7 = factor(
      age_group7,
      levels = c("0-14","15-44","45-64","65-74","75-84","85-94","95+")
    ),
    sex = factor(sex, levels = c("Female","Male"))
  )

coefs <- coef(m_qp_age7_sex)
V     <- vcov(m_qp_age7_sex)
coef_names <- names(coefs)

# ---- helpers ----

# Find coefficient name that contains exactly these parts (order-agnostic)
find_term <- function(all_names, parts) {
  parts <- as.character(parts)
  hits <- vapply(all_names, function(nm) {
    subs <- strsplit(nm, ":", fixed = TRUE)[[1]]
    all(parts %in% subs) && length(subs) == length(parts)
  }, logical(1))
  if (any(hits)) all_names[which(hits)[1]] else NA_character_
}

# Build contrast vector L (1 for selected terms, 0 otherwise)
make_L <- function(all_names, term_names) {
  L <- numeric(length(all_names)); names(L) <- all_names
  term_names <- term_names[!is.na(term_names)]
  for (tn in term_names) L[tn] <- L[tn] + 1
  L
}

base_var <- "ozone_10ppb_c"
age_fac  <- "age_group7"
sex_fac  <- "sex"

ages  <- levels(mort_env_mod[[age_fac]])  # c("0-14","15-44",..., "95+")
sexes <- levels(mort_env_mod[[sex_fac]])  # c("Female","Male")

message("Reference levels -> age_group7: ", ages[1], " ; sex: ", sexes[1])

res <- list()

for (a in ages) {
  for (s in sexes) {
    # Base ozone term
    t_base <- find_term(coef_names, c(base_var))
    
    # Two-way interactions
    t_age <- if (a != ages[1]) {
      find_term(coef_names, c(base_var, paste0(age_fac, a)))
    } else NA_character_
    
    t_sex <- if (s != sexes[1]) {
      find_term(coef_names, c(base_var, paste0(sex_fac, s)))
    } else NA_character_
    
    # Three-way interaction
    t_3way <- if (a != ages[1] && s != sexes[1]) {
      find_term(coef_names, c(base_var, paste0(age_fac, a), paste0(sex_fac, s)))
    } else NA_character_
    
    # Build contrast and compute beta, var
    L    <- make_L(coef_names, c(t_base, t_age, t_sex, t_3way))
    beta <- as.numeric(crossprod(L, coefs))
    varb <- as.numeric(t(L) %*% V %*% L)
    
    if (!is.na(varb) && varb < 0 && varb > -1e-12) varb <- 0  # numerical guard
    
    se <- sqrt(varb)
    z  <- beta / se
    p  <- 2 * (1 - pnorm(abs(z)))
    
    res[[paste(a, s, sep = "_")]] <- data.frame(
      age_group7 = a,
      sex        = s,
      logRR      = beta,
      SE         = se,
      z          = z,
      p          = p,
      RR         = exp(beta),
      LCL        = exp(beta - 1.96 * se),
      UCL        = exp(beta + 1.96 * se),
      row.names  = NULL
    )
  }
}

age_sex_rr7 <- do.call(rbind, res)
age_sex_rr7

library(ggplot2)
library(dplyr)

age_sex_rr7 <- age_sex_rr7 |>
  mutate(
    age_group7 = factor(
      age_group7,
      levels = c("0-14","15-44","45-64","65-74","75-84","85-94","95+")
    )
  )

agesex_plot <-ggplot(age_sex_rr7, aes(x = age_group7, y = RR, color = sex)) +
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  geom_errorbar(
    aes(ymin = LCL, ymax = UCL),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Age group",
    y = "Relative Risk (per 10 ppb ozone)",
    color = "Sex"
  ) +
  theme_minimal(base_size = 12)

agesex_plot

# Save to PDF
ggsave("age_mort_ox_Sex.pdf", plot = agesex_plot,
       width = 6, height = 4, units = "in")

# The age-specific risk curves show differences between male and female
# Female shows statistical significant risk decreases in age group 65-74 and 75-84, 
# whereas male shows statistical significant risk increasein age 85-94

# end.



















##### Draft #####

################################################################################
### Model with triple interaction: ozone × age × sex
m_qp_age10_sex <- glm(
  Deaths ~ ozone_10ppb_c * age_group * sex + 
    temp_c + poly(year,3) + offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)
summary(m_qp_age10_sex)

coefs <- coef(m_qp_age10_sex)
vc    <- vcov(m_qp_age10_sex)
cn    <- names(coefs)

# helper: build L contrast vector robustly (order-free matching)
find_term <- function(all, parts) {
  hits <- vapply(all, function(nm) {
    subs <- strsplit(nm, ":", fixed = TRUE)[[1]]
    all(parts %in% subs) && length(subs) == length(parts)
  }, logical(1))
  if (any(hits)) all[which(hits)[1]] else NA_character_
}

make_L <- function(all, terms) {
  L <- numeric(length(all)); names(L) <- all
  for (t in terms[!is.na(terms)]) {
    L[t] <- 1
  }
  L
}

ages  <- levels(mort_env_mod$age_group)
sexes <- levels(mort_env_mod$sex)

res <- list()
for (a in ages) {
  for (s in sexes) {
    terms <- c("ozone_10ppb_c")  # base slope
    
    if (a != ages[1]) {
      terms <- c(terms, find_term(cn, c("ozone_10ppb_c", paste0("age_group", a))))
    }
    if (s == "Male") {
      terms <- c(terms, find_term(cn, c("ozone_10ppb_c", "sexMale")))
      if (a != ages[1]) {
        terms <- c(terms, find_term(cn, c("ozone_10ppb_c", paste0("age_group", a), "sexMale")))
      }
    }
    
    L <- make_L(cn, terms)
    beta <- sum(L * coefs)
    varb <- as.numeric(t(L) %*% vc %*% L)
    if (varb < 0 && abs(varb) < 1e-12) varb <- 0  # numerical fix
    se <- sqrt(varb)
    
    z <- beta / se
    p <- 2 * pnorm(-abs(z))
    
    res[[paste(a,s,sep="_")]] <- data.frame(
      age = a, sex = s,
      logRR = beta, SE = se, z = z, p = p,
      RR = exp(beta),
      LCL = exp(beta - 1.96*se),
      UCL = exp(beta + 1.96*se)
    )
  }
}

age_sex_rr <- do.call(rbind, res)
print(age_sex_rr)
# too complex, some results likely spurious



### By age (5 age groups)
m_qp_age <- glm(
  Deaths ~ ozone_10ppb_c * age_group5 + sex + temp_c + poly(year,3) +
    offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)

summary(m_qp_age)

# --- manual RR per +10 ppb by 5 age groups (0-64 ref) ---
coefs <- coef(m_qp_age)     # <- your quasi-Poisson model object
vc    <- vcov(m_qp_age)

base_var <- "ozone_10ppb_c"     # slope name for reference age group
age_fac  <- "age_group5"        # factor name used in the model
levels5  <- c("0-64", "65-74", "75-84", "85-94", "95+")  # ensure ref first

# helpers (safe getters)
getc <- function(x, n) if (n %in% names(x)) x[[n]] else NA_real_
getv <- function(M, r, c) if (r %in% rownames(M) && c %in% colnames(M)) M[r, c] else NA_real_

# containers
k   <- length(levels5)
est <- se <- numeric(k)

# reference group ("0-64")
est[1] <- getc(coefs, base_var)
se[1]  <- sqrt(getv(vc, base_var, base_var))

# other groups: base + interaction
for (i in 2:k) {
  lev <- levels5[i]
  int_name <- paste0(base_var, ":", age_fac, lev)  # e.g., "ozone_10ppb_c:age_group565-74"
  
  beta_base <- getc(coefs, base_var)
  beta_int  <- getc(coefs, int_name)
  
  est[i] <- beta_base + beta_int
  
  var_i <- getv(vc, base_var, base_var) +
    getv(vc, int_name, int_name) +
    2 * getv(vc, base_var, int_name)
  
  se[i] <- sqrt(var_i)
}

# convert to RR and 95% CI
RR  <- exp(est)
LCL <- exp(est - 1.96 * se)
UCL <- exp(est + 1.96 * se)

# Wald z and p-values
zval <- est / se
pval <- 2 * (1 - pnorm(abs(zval)))

# tidy table
age_rr <- data.frame(
  age_group = levels5,
  logRR     = est,
  SE        = se,
  RR        = RR,
  LCL       = LCL,
  UCL       = UCL,
  z         = zval,
  p         = pval
)

age_rr

### Model with triple interaction: ozone × age × sex
m_qp_age5_sex <- glm(
  Deaths ~ ozone_10ppb_c * age_group5 * sex + 
    temp_c + year + offset(offset_log_pop),
  family = quasipoisson(link = "log"),
  data = mort_env_mod
)
summary(m_qp_age5_sex)

coefs <- coef(m_qp_age5_sex)
V     <- vcov(m_qp_age5_sex)
coef_names <- names(coefs)

# Helpers ---------------------------------------------------------------

# Find a coefficient name which contains exactly the provided 'parts' (order-agnostic)
find_term <- function(all_names, parts) {
  hits <- vapply(all_names, function(nm) {
    subs <- strsplit(nm, ":", fixed = TRUE)[[1]]
    all(parts %in% subs) && length(subs) == length(parts)
  }, logical(1))
  if (any(hits)) all_names[which(hits)[1]] else NA_character_
}

# Build a contrast vector L with 1s for the terms to sum, 0 otherwise
make_L <- function(all_names, term_names) {
  L <- numeric(length(all_names)); names(L) <- all_names
  term_names <- term_names[!is.na(term_names)]
  for (tn in term_names) L[tn] <- L[tn] + 1
  L
}

# Main variables --------------------------------------------------------
base_var <- "ozone_10ppb_c"
age_fac  <- "age_group5"
sex_fac  <- "sex"

ages  <- levels(mort_env_mod[[age_fac]])   # e.g., c("0-64","65-74","75-84","85-94","95+")
sexes <- levels(mort_env_mod[[sex_fac]])   # e.g., c("Female","Male")

# Sanity messages (which groups are references)
message("Reference levels -> age_group5: ", ages[1], " ; sex: ", sexes[1])

# Compute subgroup slopes ------------------------------------------------
res <- list()

for (a in ages) {
  for (s in sexes) {
    
    # Base ozone term
    t_base <- find_term(coef_names, c(base_var))
    
    # Two-way interactions (exist only when not at reference)
    t_age <- if (a != ages[1]) find_term(coef_names, c(base_var, paste0(age_fac, a))) else NA_character_
    t_sex <- if (s != sexes[1]) find_term(coef_names, c(base_var, paste0(sex_fac, s))) else NA_character_
    
    # Three-way interaction (only when both not reference)
    t_3way <- if (a != ages[1] && s != sexes[1]) {
      find_term(coef_names, c(base_var, paste0(age_fac, a), paste0(sex_fac, s)))
    } else NA_character_
    
    # Build contrast vector and compute beta/var
    L    <- make_L(coef_names, c(t_base, t_age, t_sex, t_3way))
    beta <- as.numeric(crossprod(L, coefs))
    varb <- as.numeric(t(L) %*% V %*% L)
    
    # guard against tiny negative round-off
    if (!is.na(varb) && varb < 0 && varb > -1e-12) varb <- 0
    
    se <- sqrt(varb)
    z  <- beta / se
    p  <- 2 * (1 - pnorm(abs(z)))
    
    res[[paste(a, s, sep = "_")]] <- data.frame(
      age_group5 = a,
      sex        = s,
      logRR      = beta,
      SE         = se,
      z          = z,
      p          = p,
      RR         = exp(beta),
      LCL        = exp(beta - 1.96 * se),
      UCL        = exp(beta + 1.96 * se),
      row.names  = NULL
    )
  }
}

age_sex_rr_5 <- do.call(rbind, res)
age_sex_rr_5
