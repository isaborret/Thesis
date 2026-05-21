# =============================================================================
# Brain-Heart Coupling: Gamma GLMM Analysis
# Author: Isa Borret
# 
# Models GC values directly with Gamma family + log link.
# Run sections sequentially; check output at each step before continuing.
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------

# Install packages
library(lme4)
library(emmeans)
library(car)
library(DHARMa)
library(MuMIn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pbkrtest)  # only needed if you run parametric bootstrap (slow)
library(patchwork)


# Set reproducibility seed
set.seed(2026) #makes any randomness reproducible, bc residuals are simulated(it generates synthetic data from your fitted model and compares them to observed) 
#With a seed, those simulations produce identical numbers each run

# Path to data file
DATA_PATH <- "C:\\Users\\Isa\\OneDrive\\Documenten\\Unif 2e master\\Masterproef_II\\data\\derivatives\\GC\\gc_results_for_R.csv"

# Whether to run parametric bootstrap LR tests (slow: ~15 min total).
# Recommend FALSE for first pass, TRUE before final write-up.
RUN_PBOOT <- TRUE
N_PBOOT   <- 1000  # bootstrap reps; 1000 is the standard


# ---- 1. Load and reshape ----------------------------------------------------

df_wide <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

cat("Subjects per stage x dream:\n")
print(with(df_wide, table(sleep_stage, dream)) / 2)  # divided by 2 because 2 electrodes per subject

# Wide -> long
gc_long <- df_wide %>%
  pivot_longer(
    cols = matches("^(EEG_to_HRV|HRV_to_EEG)_"),
    names_to = c("direction", "frequency"),
    names_pattern = "(EEG_to_HRV|HRV_to_EEG)_(.+)",
    values_to = "GC"
  )

# Add band-matched ARP_brain
gc_long <- gc_long %>%
  rowwise() %>%
  mutate(ARP_brain = case_when(
    frequency == "delta" ~ ARP_delta,
    frequency == "theta" ~ ARP_theta,
    frequency == "alpha" ~ ARP_alpha,
    frequency == "beta"  ~ ARP_beta,
    frequency == "gamma" ~ ARP_gamma
  )) %>%
  ungroup() %>%
  select(subject_id, dream, sleep_stage, electrode, frequency, direction,
         GC, ARP_brain, ARP_HRV, n_beats, morder)

# Set factors with explicit reference levels
gc_long <- gc_long %>%
  mutate(
    subject_id  = factor(subject_id),
    dream       = factor(dream, levels = c(0, 1), labels = c("no", "yes")),
    sleep_stage = factor(sleep_stage, levels = c("N2", "REM")),
    electrode   = factor(electrode, levels = c("Cz", "Pz")),
    frequency   = factor(frequency, levels = c("delta", "theta", "alpha", "beta", "gamma")),
    direction   = factor(direction, levels = c("HRV_to_EEG", "EEG_to_HRV"))
  )

# Center continuous covariates (so intercept and main effects are interpretable at the average ARP value rather than at zero)
# within each subset they'll be re-centered, but global centering
# avoids issues if you later want to compare across subsets
gc_long$ARP_brain_c <- gc_long$ARP_brain - mean(gc_long$ARP_brain, na.rm = TRUE)
gc_long$ARP_HRV_c   <- gc_long$ARP_HRV   - mean(gc_long$ARP_HRV,   na.rm = TRUE)

cat("\nLong-format dataset:\n")
print(dim(gc_long))
cat("\nSummary of GC:\n")
print(summary(gc_long$GC))
cat("\nAny zeros or non-finite:\n")
cat("  zeros:", sum(gc_long$GC == 0), "\n")
cat("  non-finite:", sum(!is.finite(gc_long$GC)), "\n")

# ---- 2. Build the 12 subsets -----------------------------------------------
# Filters out datasets per model in each approach
# Approach 1: 8 subsets (stage x electrode x direction)
# Approach 2: 4 subsets (stage x electrode)

# Approach 1 subsets
df_REM_Cz_EEGtoHRV <- gc_long %>% filter(sleep_stage=="REM", electrode=="Cz", direction=="EEG_to_HRV") %>% droplevels()
df_REM_Cz_HRVtoEEG <- gc_long %>% filter(sleep_stage=="REM", electrode=="Cz", direction=="HRV_to_EEG") %>% droplevels()
df_REM_Pz_EEGtoHRV <- gc_long %>% filter(sleep_stage=="REM", electrode=="Pz", direction=="EEG_to_HRV") %>% droplevels()
df_REM_Pz_HRVtoEEG <- gc_long %>% filter(sleep_stage=="REM", electrode=="Pz", direction=="HRV_to_EEG") %>% droplevels()
df_N2_Cz_EEGtoHRV  <- gc_long %>% filter(sleep_stage=="N2",  electrode=="Cz", direction=="EEG_to_HRV") %>% droplevels()
df_N2_Cz_HRVtoEEG  <- gc_long %>% filter(sleep_stage=="N2",  electrode=="Cz", direction=="HRV_to_EEG") %>% droplevels()
df_N2_Pz_EEGtoHRV  <- gc_long %>% filter(sleep_stage=="N2",  electrode=="Pz", direction=="EEG_to_HRV") %>% droplevels()
df_N2_Pz_HRVtoEEG  <- gc_long %>% filter(sleep_stage=="N2",  electrode=="Pz", direction=="HRV_to_EEG") %>% droplevels()

# Approach 2 subsets
df_REM_Cz <- gc_long %>% filter(sleep_stage=="REM", electrode=="Cz") %>% droplevels()
df_REM_Pz <- gc_long %>% filter(sleep_stage=="REM", electrode=="Pz") %>% droplevels()
df_N2_Cz  <- gc_long %>% filter(sleep_stage=="N2",  electrode=="Cz") %>% droplevels()
df_N2_Pz  <- gc_long %>% filter(sleep_stage=="N2",  electrode=="Pz") %>% droplevels()


# Control object reused across all glmer calls
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))


# =============================================================================
#                 APPROACH 1: SEPARATE MODELS PER DIRECTION
# Formula: GC ~ dream * frequency + (random | subject_id)
# Family: Gamma(link = "log")
# =============================================================================

cat("\n=================================================================\n")
cat("         APPROACH 1: SEPARATE MODELS PER DIRECTION               \n")
cat("=================================================================\n")


# -----------------------------------------------------------------------------
# REM | Cz | EEG -> HRV
# -----------------------------------------------------------------------------
cat("\n--- REM | Cz | EEG -> HRV | Model A: full random slope ---\n")
A1_REM_Cz_EEGtoHRV_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_REM_Cz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Cz_EEGtoHRV_full))

cat("\n--- REM | Cz | EEG -> HRV | Model B: diagonal random slope ---\n")
A1_REM_Cz_EEGtoHRV_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_REM_Cz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Cz_EEGtoHRV_diag))

cat("\n--- REM | Cz | EEG -> HRV | Model C: random intercept only ---\n")
A1_REM_Cz_EEGtoHRV_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_REM_Cz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Cz_EEGtoHRV_int))

cat("\n--- REM | Cz | EEG -> HRV | Model D: no random intercept ---\n") #this one
A1_REM_Cz_EEGtoHRV_glm <- glm(
  GC ~ 1 + dream * frequency,
  data = df_REM_Cz_EEGtoHRV, family = Gamma(link = "log")
)
print(summary(A1_REM_Cz_EEGtoHRV_glm))

# -----------------------------------------------------------------------------
# REM | Cz | HRV -> EEG
# -----------------------------------------------------------------------------
cat("\n--- REM | Cz | HRV -> EEG | Model A: full random slope ---\n")
A1_REM_Cz_HRVtoEEG_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_REM_Cz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Cz_HRVtoEEG_full))

cat("\n--- REM | Cz | HRV -> EEG | Model B: diagonal random slope ---\n")
A1_REM_Cz_HRVtoEEG_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_REM_Cz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Cz_HRVtoEEG_diag))

cat("\n--- REM | Cz | HRV -> EEG | Model C: random intercept only ---\n")
A1_REM_Cz_HRVtoEEG_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_REM_Cz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Cz_HRVtoEEG_int))

cat("\n--- REM | Cz | HRV -> EEG | Model D: no random intercept ---\n") #this one
A1_REM_Cz_HRVtoEEG_glm <- glm(
  GC ~ 1 + dream * frequency,
  data = df_REM_Cz_HRVtoEEG, family = Gamma(link = "log")
)
print(summary(A1_REM_Cz_HRVtoEEG_glm))
# -----------------------------------------------------------------------------
# REM | Pz | EEG -> HRV
# -----------------------------------------------------------------------------
cat("\n--- REM | Pz | EEG -> HRV | Model A: full random slope ---\n")
A1_REM_Pz_EEGtoHRV_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_REM_Pz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Pz_EEGtoHRV_full))

cat("\n--- REM | Pz | EEG -> HRV | Model B: diagonal random slope ---\n")
A1_REM_Pz_EEGtoHRV_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_REM_Pz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Pz_EEGtoHRV_diag))

cat("\n--- REM | Pz | EEG -> HRV | Model C: random intercept only ---\n") #this one
A1_REM_Pz_EEGtoHRV_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_REM_Pz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Pz_EEGtoHRV_int))


# -----------------------------------------------------------------------------
# REM | Pz | HRV -> EEG
# -----------------------------------------------------------------------------
cat("\n--- REM | Pz | HRV -> EEG | Model A: full random slope ---\n")
A1_REM_Pz_HRVtoEEG_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_REM_Pz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Pz_HRVtoEEG_full))

cat("\n--- REM | Pz | HRV -> EEG | Model B: diagonal random slope ---\n")
A1_REM_Pz_HRVtoEEG_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_REM_Pz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Pz_HRVtoEEG_diag))

cat("\n--- REM | Pz | HRV -> EEG | Model C: random intercept only ---\n") #this one
A1_REM_Pz_HRVtoEEG_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_REM_Pz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_REM_Pz_HRVtoEEG_int))


# -----------------------------------------------------------------------------
# N2 | Cz | EEG -> HRV
# -----------------------------------------------------------------------------
cat("\n--- N2 | Cz | EEG -> HRV | Model A: full random slope ---\n")
A1_N2_Cz_EEGtoHRV_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_N2_Cz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Cz_EEGtoHRV_full))

cat("\n--- N2 | Cz | EEG -> HRV | Model B: diagonal random slope ---\n")
A1_N2_Cz_EEGtoHRV_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_N2_Cz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Cz_EEGtoHRV_diag))

cat("\n--- N2 | Cz | EEG -> HRV | Model C: random intercept only ---\n") #this one
A1_N2_Cz_EEGtoHRV_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_N2_Cz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Cz_EEGtoHRV_int))


# -----------------------------------------------------------------------------
# N2 | Cz | HRV -> EEG
# -----------------------------------------------------------------------------
cat("\n--- N2 | Cz | HRV -> EEG | Model A: full random slope ---\n")
A1_N2_Cz_HRVtoEEG_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_N2_Cz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Cz_HRVtoEEG_full))

cat("\n--- N2 | Cz | HRV -> EEG | Model B: diagonal random slope ---\n")
A1_N2_Cz_HRVtoEEG_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_N2_Cz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Cz_HRVtoEEG_diag))

cat("\n--- N2 | Cz | HRV -> EEG | Model C: random intercept only ---\n") #this one
A1_N2_Cz_HRVtoEEG_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_N2_Cz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Cz_HRVtoEEG_int))


# -----------------------------------------------------------------------------
# N2 | Pz | EEG -> HRV
# -----------------------------------------------------------------------------
cat("\n--- N2 | Pz | EEG -> HRV | Model A: full random slope ---\n")
A1_N2_Pz_EEGtoHRV_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_N2_Pz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Pz_EEGtoHRV_full))

cat("\n--- N2 | Pz | EEG -> HRV | Model B: diagonal random slope ---\n")
A1_N2_Pz_EEGtoHRV_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_N2_Pz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Pz_EEGtoHRV_diag))

cat("\n--- N2 | Pz | EEG -> HRV | Model C: random intercept only ---\n")
A1_N2_Pz_EEGtoHRV_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_N2_Pz_EEGtoHRV, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Pz_EEGtoHRV_int))


# -----------------------------------------------------------------------------
# N2 | Pz | HRV -> EEG
# -----------------------------------------------------------------------------
cat("\n--- N2 | Pz | HRV -> EEG | Model A: full random slope ---\n")
A1_N2_Pz_HRVtoEEG_full <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency | subject_id),
  data = df_N2_Pz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Pz_HRVtoEEG_full))

cat("\n--- N2 | Pz | HRV -> EEG | Model B: diagonal random slope ---\n")
A1_N2_Pz_HRVtoEEG_diag <- glmer(
  GC ~ 1 + dream * frequency + (1 + frequency || subject_id),
  data = df_N2_Pz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Pz_HRVtoEEG_diag))

cat("\n--- N2 | Pz | HRV -> EEG | Model C: random intercept only ---\n")
A1_N2_Pz_HRVtoEEG_int <- glmer(
  GC ~ 1 + dream * frequency + (1 | subject_id),
  data = df_N2_Pz_HRVtoEEG, family = Gamma(link = "log"), control = ctrl
)
print(summary(A1_N2_Pz_HRVtoEEG_int))

cat("\n--- N2 | Pz | HRV -> EEG | Model D: no random intercept ---\n")
A1_N2_Pz_HRVtoEEG_glm <- glm(
  GC ~ 1 + dream * frequency,
  data = df_N2_Pz_HRVtoEEG, family = Gamma(link = "log")
)
print(summary(A1_N2_Pz_HRVtoEEG_glm))

# =============================================================================
#                 APPROACH 1: MODEL COMPARISONS
# =============================================================================

cat("\n=================================================================\n")
cat("         APPROACH 1: RANDOM-EFFECTS COMPARISONS                  \n")
cat("=================================================================\n")

cat("\n--- REM | Cz | EEG -> HRV ---\n")
print(anova(A1_REM_Cz_EEGtoHRV_int, A1_REM_Cz_EEGtoHRV_diag, A1_REM_Cz_EEGtoHRV_full))

cat("\n--- REM | Cz | HRV -> EEG ---\n")
print(anova(A1_REM_Cz_HRVtoEEG_int, A1_REM_Cz_HRVtoEEG_diag, A1_REM_Cz_HRVtoEEG_full))

cat("\n--- REM | Pz | EEG -> HRV ---\n")
print(anova(A1_REM_Pz_EEGtoHRV_int, A1_REM_Pz_EEGtoHRV_diag, A1_REM_Pz_EEGtoHRV_full))

cat("\n--- REM | Pz | HRV -> EEG ---\n")
print(anova(A1_REM_Pz_HRVtoEEG_int, A1_REM_Pz_HRVtoEEG_diag, A1_REM_Pz_HRVtoEEG_full))

cat("\n--- N2 | Cz | EEG -> HRV ---\n")
print(anova(A1_N2_Cz_EEGtoHRV_int, A1_N2_Cz_EEGtoHRV_diag, A1_N2_Cz_EEGtoHRV_full))

cat("\n--- N2 | Cz | HRV -> EEG ---\n")
print(anova(A1_N2_Cz_HRVtoEEG_int, A1_N2_Cz_HRVtoEEG_diag, A1_N2_Cz_HRVtoEEG_full))

cat("\n--- N2 | Pz | EEG -> HRV ---\n")
print(anova(A1_N2_Pz_EEGtoHRV_int, A1_N2_Pz_EEGtoHRV_diag, A1_N2_Pz_EEGtoHRV_full))

cat("\n--- N2 | Pz | HRV -> EEG ---\n")
print(anova(A1_N2_Pz_HRVtoEEG_int, A1_N2_Pz_HRVtoEEG_diag, A1_N2_Pz_HRVtoEEG_full))


# =============================================================================
#         APPROACH 2: COMBINED MODEL WITH DIRECTION AS PREDICTOR
# Formula: GC ~ dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
#               (random | subject_id)
# =============================================================================

cat("\n=================================================================\n")
cat("         APPROACH 2: COMBINED MODEL                              \n")
cat("=================================================================\n")


# -----------------------------------------------------------------------------
# REM | Cz
# -----------------------------------------------------------------------------
cat("\n--- REM | Cz | Model A: full random slope ---\n")
A2_REM_Cz_full <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_REM_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_REM_Cz_full))

cat("\n--- REM | Cz | Model B: diagonal random slope ---\n")
A2_REM_Cz_diag <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_REM_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_REM_Cz_diag))

cat("\n--- REM | Cz | Model C: random intercept only ---\n")
A2_REM_Cz_int <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_REM_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_REM_Cz_int))

cat("\n--- REM | Cz | Model D: no random intercept ---\n")
A2_REM_Cz_glm <- glm(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c,
  data = df_REM_Cz, family = Gamma(link = "log")
)
print(summary(A2_REM_Cz_glm))

# -----------------------------------------------------------------------------
# REM | Pz
# -----------------------------------------------------------------------------
cat("\n--- REM | Pz | Model A: full random slope ---\n")
A2_REM_Pz_full <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_REM_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_REM_Pz_full))

cat("\n--- REM | Pz | Model B: diagonal random slope ---\n")
A2_REM_Pz_diag <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_REM_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_REM_Pz_diag))

cat("\n--- REM | Pz | Model C: random intercept only ---\n")
A2_REM_Pz_int <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_REM_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_REM_Pz_int))


# -----------------------------------------------------------------------------
# N2 | Cz
# -----------------------------------------------------------------------------
cat("\n--- N2 | Cz | Model A: full random slope ---\n")
A2_N2_Cz_full <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_N2_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_N2_Cz_full))

cat("\n--- N2 | Cz | Model B: diagonal random slope ---\n")
A2_N2_Cz_diag <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_N2_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_N2_Cz_diag))

cat("\n--- N2 | Cz | Model C: random intercept only ---\n")
A2_N2_Cz_int <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_N2_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_N2_Cz_int))


# -----------------------------------------------------------------------------
# N2 | Pz
# -----------------------------------------------------------------------------
cat("\n--- N2 | Pz | Model A: full random slope ---\n")
A2_N2_Pz_full <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_N2_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_N2_Pz_full))

cat("\n--- N2 | Pz | Model B: diagonal random slope ---\n")
A2_N2_Pz_diag <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_N2_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_N2_Pz_diag))

cat("\n--- N2 | Pz | Model C: random intercept only ---\n")
A2_N2_Pz_int <- glmer(
  GC ~ 1 + dream * direction + frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_N2_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2_N2_Pz_int))


# =============================================================================
#                 APPROACH 2: MODEL COMPARISONS
# =============================================================================

cat("\n=================================================================\n")
cat("         APPROACH 2: RANDOM-EFFECTS COMPARISONS                  \n")
cat("=================================================================\n")

cat("\n--- REM | Cz ---\n")
print(anova(A2_REM_Cz_int, A2_REM_Cz_diag, A2_REM_Cz_full))

cat("\n--- REM | Pz ---\n")
print(anova(A2_REM_Pz_int, A2_REM_Pz_diag, A2_REM_Pz_full))

cat("\n--- N2 | Cz ---\n")
print(anova(A2_N2_Cz_int, A2_N2_Cz_diag, A2_N2_Cz_full))

cat("\n--- N2 | Pz ---\n")
print(anova(A2_N2_Pz_int, A2_N2_Pz_diag, A2_N2_Pz_full))


# =============================================================================
#         APPROACH 2 EXPLORATORY: THREE-WAY INTERACTION
# Formula: GC ~ dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
#               (random | subject_id)
# =============================================================================
# CHANGE NOTE: this section is now structured exactly like Approach 1 and 2,
# with all three random-effects structures fitted explicitly per subset.

cat("\n=================================================================\n")
cat("         APPROACH 2 EXPLORATORY: THREE-WAY INTERACTION           \n")
cat("=================================================================\n")


# -----------------------------------------------------------------------------
# REM | Cz | three-way
# -----------------------------------------------------------------------------
cat("\n--- REM | Cz | 3way | Model A: full random slope ---\n")
A2x_REM_Cz_full <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_REM_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_REM_Cz_full))

cat("\n--- REM | Cz | 3way | Model B: diagonal random slope ---\n")
A2x_REM_Cz_diag <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_REM_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_REM_Cz_diag))

cat("\n--- REM | Cz | 3way | Model C: random intercept only ---\n")
A2x_REM_Cz_int <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_REM_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_REM_Cz_int))

cat("\n--- REM | Cz | 3way | Model D: no random intercept ---\n")
A2x_REM_Cz_glm <- glm(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c,
  data = df_REM_Cz, family = Gamma(link = "log")
)
print(summary(A2x_REM_Cz_glm))

# -----------------------------------------------------------------------------
# REM | Pz | three-way
# -----------------------------------------------------------------------------
cat("\n--- REM | Pz | 3way | Model A: full random slope ---\n")
A2x_REM_Pz_full <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_REM_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_REM_Pz_full))

cat("\n--- REM | Pz | 3way | Model B: diagonal random slope ---\n")
A2x_REM_Pz_diag <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_REM_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_REM_Pz_diag))

cat("\n--- REM | Pz | 3way | Model C: random intercept only ---\n")
A2x_REM_Pz_int <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_REM_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_REM_Pz_int))


# -----------------------------------------------------------------------------
# N2 | Cz | three-way
# -----------------------------------------------------------------------------
cat("\n--- N2 | Cz | 3way | Model A: full random slope ---\n")
A2x_N2_Cz_full <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_N2_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_N2_Cz_full))

cat("\n--- N2 | Cz | 3way | Model B: diagonal random slope ---\n")
A2x_N2_Cz_diag <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_N2_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_N2_Cz_diag))

cat("\n--- N2 | Cz | 3way | Model C: random intercept only ---\n")
A2x_N2_Cz_int <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_N2_Cz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_N2_Cz_int))


# -----------------------------------------------------------------------------
# N2 | Pz | three-way
# -----------------------------------------------------------------------------
cat("\n--- N2 | Pz | 3way | Model A: full random slope ---\n")
A2x_N2_Pz_full <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency | subject_id),
  data = df_N2_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_N2_Pz_full))

cat("\n--- N2 | Pz | 3way | Model B: diagonal random slope ---\n")
A2x_N2_Pz_diag <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 + frequency || subject_id),
  data = df_N2_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_N2_Pz_diag))

cat("\n--- N2 | Pz | 3way | Model C: random intercept only ---\n")
A2x_N2_Pz_int <- glmer(
  GC ~ 1 + dream * direction * frequency + ARP_brain_c + ARP_HRV_c +
    (1 | subject_id),
  data = df_N2_Pz, family = Gamma(link = "log"), control = ctrl
)
print(summary(A2x_N2_Pz_int))


# =============================================================================
#                 APPROACH 2 EXPLORATORY: MODEL COMPARISONS
# =============================================================================

cat("\n=================================================================\n")
cat("         APPROACH 2 EXPLORATORY: COMPARISONS                     \n")
cat("=================================================================\n")

# Random-effects comparisons within the three-way structure
cat("\n--- REM | Cz | 3way: RE comparison ---\n")
print(anova(A2x_REM_Cz_int, A2x_REM_Cz_diag, A2x_REM_Cz_full))

cat("\n--- REM | Pz | 3way: RE comparison ---\n")
print(anova(A2x_REM_Pz_int, A2x_REM_Pz_diag, A2x_REM_Pz_full))

cat("\n--- N2 | Cz | 3way: RE comparison ---\n")
print(anova(A2x_N2_Cz_int, A2x_N2_Cz_diag, A2x_N2_Cz_full))

cat("\n--- N2 | Pz | 3way: RE comparison ---\n")
print(anova(A2x_N2_Pz_int, A2x_N2_Pz_diag, A2x_N2_Pz_full))

# Additive (Approach 2 main) vs three-way (Approach 2 exploratory),
# at random-intercept-only structure for fair comparison
cat("\n--- REM | Cz | additive vs three-way ---\n")
print(anova(A2_REM_Cz_glm, A2x_REM_Cz_glm))

cat("\n--- REM | Pz | additive vs three-way ---\n")
print(anova(A2_REM_Pz_int, A2x_REM_Pz_int))

cat("\n--- N2 | Cz | additive vs three-way ---\n")
print(anova(A2_N2_Cz_int, A2x_N2_Cz_int))

cat("\n--- N2 | Pz | additive vs three-way ---\n")
print(anova(A2_N2_Pz_int, A2x_N2_Pz_int))


# =============================================================================
# Once you've inspected the model summaries and comparisons, MANUALLY assign
# the best-fitting model from each subset to the variables below. The rest of
# the script (inference, effect sizes, diagnostics) uses these.
#
# Selection rule: most complex non-singular model with lowest AIC.
# If a model is singular (variance estimated at 0) or fails to converge,
# don't pick it. Default to _int if neither richer model fits cleanly.
# =============================================================================

cat("\n=================================================================\n")
cat("         BEST-MODEL ASSIGNMENT (edit manually)                   \n")
cat("=================================================================\n")

# Approach 1
A1_REM_Cz_EEGtoHRV_best <- A1_REM_Cz_EEGtoHRV_glm   # <-- edit
A1_REM_Cz_HRVtoEEG_best <- A1_REM_Cz_HRVtoEEG_glm   # <-- edit
A1_REM_Pz_EEGtoHRV_best <- A1_REM_Pz_EEGtoHRV_int   # <-- edit
A1_REM_Pz_HRVtoEEG_best <- A1_REM_Pz_HRVtoEEG_int   # <-- edit
A1_N2_Cz_EEGtoHRV_best  <- A1_N2_Cz_EEGtoHRV_int    # <-- edit
A1_N2_Cz_HRVtoEEG_best  <- A1_N2_Cz_HRVtoEEG_int    # <-- edit
A1_N2_Pz_EEGtoHRV_best  <- A1_N2_Pz_EEGtoHRV_int    # <-- edit
A1_N2_Pz_HRVtoEEG_best  <- A1_N2_Pz_HRVtoEEG_glm    # <-- edit

# Approach 2
A2_REM_Cz_best <- A2_REM_Cz_glm   # <-- edit
A2_REM_Pz_best <- A2_REM_Pz_int   # <-- edit
A2_N2_Cz_best  <- A2_N2_Cz_int    # <-- edit
A2_N2_Pz_best  <- A2_N2_Pz_int    # <-- edit

# Approach 2 exploratory
A2x_REM_Cz_best <- A2x_REM_Cz_glm   # <-- edit
A2x_REM_Pz_best <- A2x_REM_Pz_int   # <-- edit
A2x_N2_Cz_best  <- A2x_N2_Cz_int    # <-- edit
A2x_N2_Pz_best  <- A2x_N2_Pz_int    # <-- edit

# =============================================================================
#                 PARAMETRIC BOOTSTRAP / LR TESTS — ALL FIXED-EFFECT TERMS
# =============================================================================
# For each glmer model: PBmodcomp on every fixed-effect term.
# For each glm fallback: anova(test="Chisq") on every fixed-effect term.
#
# Marginality principle: when testing a main effect, drop both the main
# effect AND every interaction it appears in.
#
# Storage: each test stores its p-value in `pb_results` with a label like
# "A1_REM_Cz_EEGtoHRV | dream" so the FDR correction at the end can pick
# out the dream-related ones.
# =============================================================================

# Helper: run a test and return its p-value.
# Detects whether the model is glmer (use PBmodcomp) or glm (use anova LR).
test_term <- function(model_full, model_red, label, n_pboot = N_PBOOT) {
  # Skip if already done
  if (!is.null(pb_results[[label]])) {
    cat("\n--- Skipping", label, "(already done) ---\n")
    return(pb_results[[label]])
  }
  cat("\n---", label, "---\n")
  if (inherits(model_full, "merMod")) {
    res <- suppressWarnings(PBmodcomp(model_full, model_red, nsim = n_pboot))
    print(res)
    p <- res$test$p.value[2]
  } else {
    res <- anova(model_full, model_red, test = "Chisq")
    print(res)
    p <- res$`Pr(>Chi)`[2]
  }
  pb_results[[label]] <<- p
  save(pb_results, file = "pb_results_running.RData")
  p
}

if (!exists("pb_results")) pb_results <- list()


if (RUN_PBOOT) {
  
  cat("\n=================================================================\n")
  cat("         PARAMETRIC BOOTSTRAP / LR TESTS                          \n")
  cat("=================================================================\n")
  
  # ---- Approach 1: each subset has terms dream, frequency, dream:frequency ----
  
  for (key in c("A1_REM_Cz_EEGtoHRV", "A1_REM_Cz_HRVtoEEG",
                "A1_REM_Pz_EEGtoHRV", "A1_REM_Pz_HRVtoEEG",
                "A1_N2_Cz_EEGtoHRV",  "A1_N2_Cz_HRVtoEEG",
                "A1_N2_Pz_EEGtoHRV",  "A1_N2_Pz_HRVtoEEG")) {
    
    m <- get(paste0(key, "_best"))
    
    # dream:frequency interaction
    pb_results[[paste(key, "dream:frequency", sep=" | ")]] <-
      test_term(m, update(m, . ~ . - dream:frequency),
                paste(key, "dream:frequency"))
    
    # dream main (drop dream + dream:frequency)
    pb_results[[paste(key, "dream", sep=" | ")]] <-
      test_term(m, update(m, . ~ . - dream - dream:frequency),
                paste(key, "dream"))
    
    # frequency main (drop frequency + dream:frequency)
    pb_results[[paste(key, "frequency", sep=" | ")]] <-
      test_term(m, update(m, . ~ . - frequency - dream:frequency),
                paste(key, "frequency"))
  }
  
  save(pb_results, file = "pb_results_after_A1.RData")
  cat("\n>> Approach 1 PB tests complete. Saved to pb_results_after_A1.RData\n")
  
  
  # ---- Approach 2: terms dream, direction, dream:direction, frequency,
  #                  ARP_brain_c, ARP_HRV_c ----
  
  for (key in c("A2_REM_Cz", "A2_REM_Pz", "A2_N2_Cz", "A2_N2_Pz")) {
    
    m <- get(paste0(key, "_best"))
    
    # dream:direction interaction
    test_term(m, update(m, . ~ . - dream:direction),
                paste(key, "dream:direction", sep=" | "))
    
    # dream main (drop dream + dream:direction)
    test_term(m, update(m, . ~ . - dream - dream:direction),
                paste(key, "dream", sep=" | "))
    
    # direction main (drop direction + dream:direction)
    test_term(m, update(m, . ~ . - direction - dream:direction),
                paste(key, "direction", sep=" | "))
    
    # frequency main
    test_term(m, update(m, . ~ . - frequency),
                paste(key, "frequency", sep=" | "))
    
    # ARP_brain_c
    test_term(m, update(m, . ~ . - ARP_brain_c),
                paste(key, "ARP_brain_c", sep=" | "))
    
    # ARP_HRV_c
    test_term(m, update(m, . ~ . - ARP_HRV_c),
                paste(key, "ARP_HRV_c", sep=" | "))
  }
  
  save(pb_results, file = "pb_results_after_A2.RData")
  cat("\n>> Approach 2 PB tests complete. Saved to pb_results_after_A2.RData\n")
  
  
  # ---- Approach 2 exploratory: terms with three-way structure ----
  # All factors interact: dream * direction * frequency
  # Tests of main effects must drop all interactions involving that term.
  
  for (key in c("A2x_REM_Cz", "A2x_REM_Pz", "A2x_N2_Cz", "A2x_N2_Pz")) {
    
    m <- get(paste0(key, "_best"))
    
    # three-way interaction
    test_term(m, update(m, . ~ . - dream:direction:frequency),
                paste(key, "dream:direction:frequency", sep=" | "))
    
    # dream:direction (drop dream:direction + the three-way)
    test_term(m, update(m, . ~ . - dream:direction - dream:direction:frequency),
                paste(key, "dream:direction", sep=" | "))
    
    # dream:frequency (drop dream:frequency + the three-way)
    test_term(m, update(m, . ~ . - dream:frequency - dream:direction:frequency),
                paste(key, "dream:frequency", sep=" | "))
    
    # direction:frequency (drop direction:frequency + the three-way)
    test_term(m, update(m, . ~ . - direction:frequency - dream:direction:frequency),
                paste(key, "direction:frequency", sep=" | "))
    
    # dream main (drop dream + dream:direction + dream:frequency + 3-way)
    test_term(m, update(m, . ~ . - dream - dream:direction - dream:frequency - dream:direction:frequency),
                paste(key, "dream", sep=" | "))
    
    # direction main
    test_term(m, update(m, . ~ . - direction - dream:direction - direction:frequency - dream:direction:frequency),
                paste(key, "direction", sep=" | "))
    
    # frequency main
    test_term(m, update(m, . ~ . - frequency - dream:frequency - direction:frequency - dream:direction:frequency),
                paste(key, "frequency", sep=" | "))
    
    # ARP_brain_c
    test_term(m, update(m, . ~ . - ARP_brain_c),
                paste(key, "ARP_brain_c", sep=" | "))
    
    # ARP_HRV_c
    test_term(m, update(m, . ~ . - ARP_HRV_c),
                paste(key, "ARP_HRV_c", sep=" | "))
  }
  
  save(pb_results, file = "pb_results_complete.RData")
  cat("\n>> All PB tests complete. Saved to pb_results_complete.RData\n")
}


# =============================================================================
#                 RESULTS TABLE: all PB tests, sorted
# =============================================================================

pb_table <- data.frame(
  contrast = names(pb_results),
  p_raw    = unlist(pb_results),
  stringsAsFactors = FALSE
)
pb_table <- pb_table[order(pb_table$p_raw), ]

cat("\n=================================================================\n")
cat("         ALL PARAMETRIC BOOTSTRAP / LR P-VALUES                    \n")
cat("=================================================================\n")
print(pb_table, row.names = FALSE)


# =============================================================================
#                 FDR CORRECTION on all dream-related tests
# =============================================================================
# Identifying dream-related contrasts:
# - any term containing "dream" (dream main, dream:frequency, dream:direction,
#   dream:direction:frequency)

is_dream_related <- grepl("dream", pb_table$contrast)
focal_table <- pb_table[is_dream_related, ]
focal_table$p_FDR <- p.adjust(focal_table$p_raw, method = "BH")

cat("\n=================================================================\n")
cat("         FOCAL TESTS (DREAM-RELATED) WITH FDR CORRECTION           \n")
cat("=================================================================\n")
cat("Number of focal tests:", nrow(focal_table), "\n\n")
print(focal_table[order(focal_table$p_FDR), ], row.names = FALSE)


# Non-focal: reported but not FDR-corrected
nonfocal_table <- pb_table[!is_dream_related, ]
cat("\n--- Non-focal tests (descriptive, not FDR-corrected) ---\n")
print(nonfocal_table, row.names = FALSE)


cat("\n=================================================================\n")
cat("         BUILD RESULTS TABLE           \n")
cat("=================================================================\n")

# Build the results table
results <- data.frame(
  contrast = names(pb_results),
  p_raw = unlist(pb_results),
  stringsAsFactors = FALSE
)

# FDR for dream-related tests only
results$p_FDR <- NA
focal_idx <- grepl("dream", results$contrast)
results$p_FDR[focal_idx] <- p.adjust(results$p_raw[focal_idx], method = "BH")

# Parse model and term from the contrast label
results$model <- sapply(results$contrast, function(label) {
  trimws(strsplit(label, "\\|")[[1]][1])
})

results$term <- sapply(results$contrast, function(label) {
  parts <- strsplit(label, "\\|")[[1]]
  term <- trimws(parts[2])
  gsub(" \\(Wald\\)$", "", term)  # strip the (Wald) suffix
})

# Define a fixed term order so each model's rows are consistent
term_order <- c("dream", "direction", "frequency",
                "dream:direction", "dream:frequency", "direction:frequency",
                "dream:direction:frequency",
                "ARP_brain_c", "ARP_HRV_c")

# Define a model order
model_order <- c(
  "A1_REM_Cz_EEGtoHRV", "A1_REM_Cz_HRVtoEEG",
  "A1_REM_Pz_EEGtoHRV", "A1_REM_Pz_HRVtoEEG",
  "A1_N2_Cz_EEGtoHRV",  "A1_N2_Cz_HRVtoEEG",
  "A1_N2_Pz_EEGtoHRV",  "A1_N2_Pz_HRVtoEEG",
  "A2_REM_Cz", "A2_REM_Pz", "A2_N2_Cz", "A2_N2_Pz",
  "A2x_REM_Cz", "A2x_REM_Pz", "A2x_N2_Cz", "A2x_N2_Pz"
)

results$model <- factor(results$model, levels = model_order)
results$term  <- factor(results$term,  levels = term_order)

# Sort: model first, then term within model
results <- results[order(results$model, results$term), ]

# Round p-values for readability
results$p_raw <- round(results$p_raw, 4)
results$p_FDR <- round(results$p_FDR, 4)

# Reorder columns: model, term, p_raw, p_FDR
results <- results[, c("model", "term", "p_raw", "p_FDR")]

# Save
write.csv(results, "results_summary2.csv", row.names = FALSE)




# =============================================================================
#                 EFFECT SIZES: rate ratios via emmeans
# (unchanged from previous version — see gc_glmm_full.R)
# =============================================================================

# Approach 1: dream contrast per band, plus marginal dream contrast
# Approach 2: dream x direction interaction contrasts
# Approach 2 exploratory: dream x direction x frequency contrasts

# ============================================================================
# 1. THE FDR-SIGNIFICANT EFFECT
# ============================================================================
# A1_N2_Pz_HRVtoEEG: dream main effect (the only FDR-significant finding)
# Model: GC ~ dream * frequency  (this is a glm)

m <- A1_N2_Pz_HRVtoEEG_best

# Marginal dream contrast averaged across frequency bands
emm <- emmeans(m, ~ dream, type = "response")
cat("\n--- A1_N2_Pz_HRVtoEEG: marginal means by dream group ---\n")
print(emm)

cat("\n--- A1_N2_Pz_HRVtoEEG: dream contrast (yes vs no) ---\n")
print(pairs(emm, reverse = TRUE, infer = c(TRUE, TRUE)))

# Per-band breakdown — useful for understanding which bands drive the effect
emm_band <- emmeans(m, ~ dream | frequency, type = "response")
cat("\n--- A1_N2_Pz_HRVtoEEG: dream contrast per frequency band ---\n")
print(pairs(emm_band, reverse = TRUE, infer = c(TRUE, TRUE)))


# ============================================================================
# 2. THE MARGINAL FOCAL EFFECTS (raw-significant, not FDR-significant)
# ============================================================================

# A2_N2_Pz: dream main effect AND dream:direction interaction
# Model: GC ~ dream * direction + frequency + ARP_brain_c + ARP_HRV_c (glmer)

m <- A2_N2_Pz_best

# Marginal means for the 4 dream × direction cells
emm <- emmeans(m, ~ dream * direction, type = "response")
cat("\n--- A2_N2_Pz: marginal means by dream × direction ---\n")
print(emm)

# Hypothesis 1 + 2 contrast: dream effect within each direction
cat("\n--- A2_N2_Pz: dream contrast within each direction ---\n")
print(pairs(emm, by = "direction", reverse = TRUE, type = "response"))

# Hypothesis 3 contrast: direction effect within each dream group
cat("\n--- A2_N2_Pz: direction contrast within each dream group ---\n")
print(pairs(emm, by = "dream", type = "response"))

# Marginal dream contrast averaged across directions (the dream main effect)
emm_dream <- emmeans(m, ~ dream, type = "response")
cat("\n--- A2_N2_Pz: dream main effect (averaged across directions) ---\n")
print(pairs(emm_dream, reverse = TRUE, infer = c(TRUE, TRUE)))

# =============================================================================
#                 EMMEANS FOR EVERY
# =============================================================================

get_marginal_contrasts <- function(model, factors, label) {
  out <- list()
  for (f in factors) {
    if (f %in% all.vars(formula(model))) {
      formula_str <- as.formula(paste("~", f))
      emm <- emmeans(model, formula_str, type = "response")
      ct <- pairs(emm, reverse = TRUE, infer = c(TRUE, TRUE))
      ct_df <- as.data.frame(ct)
      ct_df$model <- label
      ct_df$factor <- f
      out[[f]] <- ct_df
    }
  }
  out
}

# Example for one model
result <- get_marginal_contrasts(A2_N2_Cz_best, 
                                 c("dream", "direction", "frequency"),
                                 "A2_N2_Cz")
print(result$dream)
print(result$direction)
print(result$frequency)
# =============================================================================
#                 NAKAGAWA R^2 (glmer models only)
# =============================================================================
# For glm fallbacks, use 1 - deviance/null.deviance as a pseudo-R2.

#NAKAGAWA's R² for model fit

r2_results <- data.frame()

# glmer models — Nakagawa
glmer_keys <- c("A1_REM_Pz_EEGtoHRV", "A1_REM_Pz_HRVtoEEG",
                "A1_N2_Cz_EEGtoHRV",  "A1_N2_Cz_HRVtoEEG",  "A1_N2_Pz_EEGtoHRV",
                "A2_REM_Pz", "A2_N2_Cz", "A2_N2_Pz",
                "A2x_REM_Pz", "A2x_N2_Cz", "A2x_N2_Pz")

for (key in paste0(glmer_keys, "_best")) {
  m <- get(key)
  r2 <- tryCatch(r.squaredGLMM(m), error = function(e) NULL)
  if (!is.null(r2)) {
    r2_results <- rbind(r2_results, data.frame(
      model = key,
      type = "Nakagawa",
      R2_marginal = r2[1, "R2m"],
      R2_conditional = r2[1, "R2c"]
    ))
  } else {
    r2_results <- rbind(r2_results, data.frame(
      model = key, type = "Nakagawa",
      R2_marginal = NA, R2_conditional = NA
    ))
  }
}

# glm models — pseudo-R² (1 - residual deviance / null deviance)
glm_keys <- c("A1_REM_Cz_EEGtoHRV", "A1_REM_Cz_HRVtoEEG", "A1_N2_Pz_HRVtoEEG",
              "A2_REM_Cz", "A2x_REM_Cz")

for (key in paste0(glm_keys, "_best")) {
  m <- get(key)
  pseudo_r2 <- 1 - m$deviance / m$null.deviance
  r2_results <- rbind(r2_results, data.frame(
    model = key, type = "pseudo-R2",
    R2_marginal = pseudo_r2, R2_conditional = NA
  ))
}

print(r2_results)

write.csv(r2_results, "r2_summary.csv", row.names = FALSE)

# ============================================================================
# PLOTTING OF EMMEANS
# ============================================================================



# ============================================================================
# PLOTTING OF RATE RATIOS
# ============================================================================

# ---- Helper: extract dream rate ratio + CI ---------------------------------
get_dream_rr <- function(model, by = NULL, label) {
  emm <- emmeans(model, ~ dream, by = by, type = "response")
  ct <- as.data.frame(pairs(emm, reverse = TRUE, infer = c(TRUE, TRUE)))
  
  # Handle different CI column names (glmer = asymp.LCL, glm = lower.CL)
  lower_col <- if ("asymp.LCL" %in% names(ct)) "asymp.LCL" else "lower.CL"
  upper_col <- if ("asymp.UCL" %in% names(ct)) "asymp.UCL" else "upper.CL"
  
  out <- data.frame(
    label = label,
    ratio = ct$ratio,
    lower = ct[[lower_col]],
    upper = ct[[upper_col]],
    p     = ct$p.value
  )
  
  # Add the "by" column if present (for A2 plot)
  if (!is.null(by) && by %in% names(ct)) {
    out$direction <- as.character(ct[[by]])
  }
  
  out
}


# ============================================================================
# PLOT 1: APPROACH 1 — DREAM RATE RATIOS ACROSS 8 SUBSETS
# ============================================================================

a1_rr <- rbind(
  get_dream_rr(A1_REM_Cz_EEGtoHRV_best, label = "REM | Cz | Brain→Heart"),
  get_dream_rr(A1_REM_Cz_HRVtoEEG_best, label = "REM | Cz | Heart→Brain"),
  get_dream_rr(A1_REM_Pz_EEGtoHRV_best, label = "REM | Pz | Brain→Heart"),
  get_dream_rr(A1_REM_Pz_HRVtoEEG_best, label = "REM | Pz | Heart→Brain"),
  get_dream_rr(A1_N2_Cz_EEGtoHRV_best,  label = "N2 | Cz | Brain→Heart"),
  get_dream_rr(A1_N2_Cz_HRVtoEEG_best,  label = "N2 | Cz | Heart→Brain"),
  get_dream_rr(A1_N2_Pz_EEGtoHRV_best,  label = "N2 | Pz | Brain→Heart"),
  get_dream_rr(A1_N2_Pz_HRVtoEEG_best,  label = "N2 | Pz | Heart→Brain")
)

# Add stage and direction columns for facetting/colouring
a1_rr$stage <- ifelse(grepl("REM", a1_rr$label), "REM", "N2")
a1_rr$direction <- ifelse(grepl("Brain→Heart", a1_rr$label), 
                          "Brain→Heart", "Heart→Brain")
a1_rr$stage <- factor(a1_rr$stage, levels = c("REM", "N2"))
a1_rr$label <- factor(a1_rr$label, levels = rev(a1_rr$label))  # top-down order

p_a1 <- ggplot(a1_rr, aes(x = ratio, y = label, xmin = lower, xmax = upper,
                          color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(height = 0.2) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("Brain→Heart" = "#D55E00", 
                                "Heart→Brain" = "#0072B2")) +
  labs(x = "Rate ratio: dream / no-dream (log scale)",
       y = NULL,
       title = "Approach 1: dream effect on GC across all subsets",
       color = "Direction") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor.x = element_blank())

print(p_a1)
ggsave("plot_A1_dream_rate_ratios.png", p_a1, width = 8, height = 5, dpi = 300)


# ============================================================================
# PLOT 2: APPROACH 2 — DREAM RATE RATIOS WITHIN EACH DIRECTION
# ============================================================================

a2_rr <- rbind(
  get_dream_rr(A2_REM_Cz_best, by = "direction", label = "REM | Cz"),
  get_dream_rr(A2_REM_Pz_best, by = "direction", label = "REM | Pz"),
  get_dream_rr(A2_N2_Cz_best,  by = "direction", label = "N2 | Cz"),
  get_dream_rr(A2_N2_Pz_best,  by = "direction", label = "N2 | Pz")
)

# Recode direction with arrows for consistency with Plot 1
a2_rr$direction <- ifelse(a2_rr$direction == "EEG_to_HRV",
                          "Brain→Heart", "Heart→Brain")
a2_rr$stage <- ifelse(grepl("REM", a2_rr$label), "REM", "N2")
a2_rr$stage <- factor(a2_rr$stage, levels = c("REM", "N2"))
a2_rr$label <- factor(a2_rr$label, levels = rev(unique(a2_rr$label)))

p_a2 <- ggplot(a2_rr, aes(x = ratio, y = label, xmin = lower, xmax = upper,
                          color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(height = 0.2, position = position_dodge(width = 0.4)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_x_log10() +
  scale_color_manual(values = c("Brain→Heart" = "#D55E00", 
                                "Heart→Brain" = "#0072B2")) +
  labs(x = "Rate ratio: dream / no-dream (log scale)",
       y = NULL,
       title = "Approach 2: dream effect within each direction",
       color = "Direction") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor.x = element_blank())

print(p_a2)
ggsave("plot_A2_dream_rate_ratios.png", p_a2, width = 8, height = 4, dpi = 300)

# =============================================================================
#                 MAKE COEFFICIENTS TABLE
# =============================================================================

# Extract dream-related coefficients from all models
extract_coefs <- function(model, label) {
  cf <- summary(model)$coefficients
  
  # Find rows containing "dream"
  dream_rows <- cf[grepl("dream", rownames(cf)), , drop = FALSE]
  
  if (nrow(dream_rows) == 0) return(NULL)
  
  data.frame(
    model = label,
    term = rownames(dream_rows),
    beta = round(dream_rows[, 1], 3),
    SE = round(dream_rows[, 2], 3),
    z_or_t = round(dream_rows[, 3], 3),
    p_wald = round(dream_rows[, 4], 4),
    stringsAsFactors = FALSE
  )
}

# Apply across all models
all_coefs <- do.call(rbind, list(
  extract_coefs(A1_REM_Cz_EEGtoHRV_best, "A1_REM_Cz_EEGtoHRV"),
  extract_coefs(A1_REM_Cz_HRVtoEEG_best, "A1_REM_Cz_HRVtoEEG"),
  extract_coefs(A1_REM_Pz_EEGtoHRV_best, "A1_REM_Pz_EEGtoHRV"),
  extract_coefs(A1_REM_Pz_HRVtoEEG_best, "A1_REM_Pz_HRVtoEEG"),
  extract_coefs(A1_N2_Cz_EEGtoHRV_best,  "A1_N2_Cz_EEGtoHRV"),
  extract_coefs(A1_N2_Cz_HRVtoEEG_best,  "A1_N2_Cz_HRVtoEEG"),
  extract_coefs(A1_N2_Pz_EEGtoHRV_best,  "A1_N2_Pz_EEGtoHRV"),
  extract_coefs(A1_N2_Pz_HRVtoEEG_best,  "A1_N2_Pz_HRVtoEEG"),
  extract_coefs(A2_REM_Cz_best, "A2_REM_Cz"),
  extract_coefs(A2_REM_Pz_best, "A2_REM_Pz"),
  extract_coefs(A2_N2_Cz_best,  "A2_N2_Cz"),
  extract_coefs(A2_N2_Pz_best,  "A2_N2_Pz"),
  extract_coefs(A2x_REM_Cz_best, "A2x_REM_Cz"),
  extract_coefs(A2x_REM_Pz_best, "A2x_REM_Pz"),
  extract_coefs(A2x_N2_Cz_best,  "A2x_N2_Cz"),
  extract_coefs(A2x_N2_Pz_best,  "A2x_N2_Pz")
))

write.csv(all_coefs, "dream_coefficients_all_models.csv", row.names = FALSE)
# =============================================================================
#                 MOST RECENT PLOTSSSSSSSSSSSSSS
# =============================================================================
# =============================================================================
#                 REVISED PLOTS
#   - Raw data points shown (jittered) instead of CI error bars
#   - Emmeans predicted means overlaid as larger summary points
#   - No connecting lines (x-axis categories are nominal)
# =============================================================================

# ── helpers ──────────────────────────────────────────────────────────────────

# Extract raw data from a model
# Gamma GLM with log link: model.frame() returns the response on the original
# scale already — the log link is applied internally, so no back-transform needed
get_raw_df <- function(model) {
  df <- model.frame(model)
  resp_col <- names(df)[1]
  df$GC_raw <- df[[resp_col]]   # already on original scale
  df
}

# Extract emmeans summary (means only, no CIs)
get_emm_df <- function(model, specs) {
  emm <- emmeans(model, specs = specs, type = "response")
  as.data.frame(emm)
}

# Common theme
theme_panel <- theme_bw(base_size = 9) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank()
  )

dream_colors <- c("no" = "#FF1493", "yes" = "#008000")   # coral/pink vs teal/green
dream_labels  <- c("No recall", "Dream recall")
freq_levels   <- c("delta", "theta", "alpha", "beta", "gamma")


# =============================================================================
# FIGURE 1 — APPROACH 1
# 8 panels (stage × electrode × direction)
# x = frequency band, colour = dream group
# Raw points jittered; emmeans mean overlaid as larger filled point
# =============================================================================

plot_a1 <- function(model, title) {
  
  raw <- get_raw_df(model)
  raw$frequency <- factor(raw$frequency, levels = freq_levels)
  
  emm <- get_emm_df(model, ~ dream * frequency)
  emm$frequency <- factor(emm$frequency, levels = freq_levels)
  
  ggplot() +
    # ── raw data ──────────────────────────────────────────────────────────
    geom_jitter(
      data   = raw,
      aes(x = frequency, y = GC_raw, colour = dream),
      width  = 0.18, height = 0, alpha = 0.6, size = 1.2,
      show.legend = FALSE
    ) +
    # ── emmeans means ─────────────────────────────────────────────────────
    geom_point(
      data     = emm,
      aes(x = frequency, y = response, colour = dream),
      size     = 3.5, shape = 18,          # diamond = "summary" symbol
      position = position_dodge(0.4)
    ) +
    scale_colour_manual(values = dream_colors, labels = dream_labels,
                        name = NULL) +
    scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.04)) +
    labs(x = "Frequency band", y = "GC", title = title) +
    theme_panel
}

p1_1 <- plot_a1(A1_REM_Cz_EEGtoHRV_best, "REM | Cz | Brain\u2192Heart")
p1_2 <- plot_a1(A1_REM_Cz_HRVtoEEG_best, "REM | Cz | Heart\u2192Brain")
p1_3 <- plot_a1(A1_REM_Pz_EEGtoHRV_best, "REM | Pz | Brain\u2192Heart")
p1_4 <- plot_a1(A1_REM_Pz_HRVtoEEG_best, "REM | Pz | Heart\u2192Brain")
p1_5 <- plot_a1(A1_N2_Cz_EEGtoHRV_best,  "N2  | Cz | Brain\u2192Heart")
p1_6 <- plot_a1(A1_N2_Cz_HRVtoEEG_best,  "N2  | Cz | Heart\u2192Brain")
p1_7 <- plot_a1(A1_N2_Pz_EEGtoHRV_best,  "N2  | Pz | Brain\u2192Heart")
p1_8 <- plot_a1(A1_N2_Pz_HRVtoEEG_best,  "N2  | Pz | Heart\u2192Brain")

fig1 <- (p1_1 | p1_2 | p1_3 | p1_4) /
  (p1_5 | p1_6 | p1_7 | p1_8) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title  = "Approach 1: GC by frequency band and dream recall",
    caption = "Points: individual observations (jittered). Diamonds: group means (emmeans).",
    theme  = theme(
      plot.title   = element_text(face = "bold"),
      plot.caption = element_text(size = 7),
      legend.position = "bottom"
    )
  )

print(fig1)
ggsave("fig1_approach1_scatter.png", fig1, width = 14, height = 6, dpi = 300)


# =============================================================================
# FIGURE 2 — APPROACH 2
# 4 panels (stage × electrode)
# x = direction (Brain→Heart vs Heart→Brain), colour = dream group
# =============================================================================

relabel_dir <- function(x) {
  ifelse(x == "EEG_to_HRV", "Brain\u2192Heart", "Heart\u2192Brain")
}

plot_a2 <- function(model, title) {
  
  raw <- get_raw_df(model)
  raw$direction <- factor(relabel_dir(raw$direction),
                          levels = c("Brain\u2192Heart", "Heart\u2192Brain"))
  
  emm <- get_emm_df(model, ~ dream * direction)
  emm$direction <- factor(relabel_dir(emm$direction),
                          levels = c("Brain\u2192Heart", "Heart\u2192Brain"))
  
  ggplot() +
    geom_jitter(
      data   = raw,
      aes(x = direction, y = GC_raw, colour = dream),
      width  = 0.15, height = 0, alpha = 0.6, size = 1.2,
      show.legend = FALSE
    ) +
    geom_point(
      data     = emm,
      aes(x = direction, y = response, colour = dream),
      size     = 3.5, shape = 18,
      position = position_dodge(0.4)
    ) +
    scale_colour_manual(values = dream_colors, labels = dream_labels,
                        name = NULL) +
    scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.04)) +
    labs(x = "Direction", y = "GC", title = title) +
    theme_panel
}

p2_1 <- plot_a2(A2_REM_Cz_best, "REM | Cz")
p2_2 <- plot_a2(A2_REM_Pz_best, "REM | Pz")
p2_3 <- plot_a2(A2_N2_Cz_best,  "N2  | Cz")
p2_4 <- plot_a2(A2_N2_Pz_best,  "N2  | Pz")

fig2 <- (p2_1 | p2_2 | p2_3 | p2_4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = "Approach 2: GC by direction and dream recall",
    caption = "Points: individual observations (jittered). Diamonds: group means (emmeans).",
    theme   = theme(
      plot.title   = element_text(face = "bold"),
      plot.caption = element_text(size = 7),
      legend.position = "bottom"
    )
  )

print(fig2)
ggsave("fig2_approach2_scatter.png", fig2, width = 10, height = 4, dpi = 300)

# FIGURE 2 — APPROACH 2 TO MIMIC PANELS 
# 4 panels (stage × electrode)
# x = dream recall group, colour + linetype = direction
# Brain→Heart: solid blue line; Heart→Brain: dashed orange line
# Raw points jittered; emmeans means connected by lines
# =============================================================================
relabel_dir <- function(x) {
  ifelse(x == "EEG_to_HRV", "Brain\u2192Heart", "Heart\u2192Brain")
}

dir_colors   <- c("Brain\u2192Heart" = "#0072B2", "Heart\u2192Brain" = "#ff7000")
dir_linetypes <- c("Brain\u2192Heart" = "solid",   "Heart\u2192Brain" = "dashed")
dream_x_labels <- c("yes" = "Dream recall", "no" = "No recall")

plot_a2 <- function(model, title) {
  
  raw <- get_raw_df(model)
  raw$direction <- factor(relabel_dir(raw$direction),
                          levels = c("Brain\u2192Heart", "Heart\u2192Brain"))
  raw$dream <- factor(raw$dream, levels = c("yes", "no"))
  
  emm <- get_emm_df(model, ~ dream * direction)
  emm$direction <- factor(relabel_dir(emm$direction),
                          levels = c("Brain\u2192Heart", "Heart\u2192Brain"))
  emm$dream <- factor(emm$dream, levels = c("yes", "no"))
  
  ggplot() +
    # raw data: jittered points, dodged by direction
    geom_point(
      data     = raw,
      aes(x = dream, y = GC_raw, colour = direction),
      alpha    = 0.5, size = 1.2,
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.4),
      show.legend = FALSE
    ) +
    # emmeans means connected by lines
    geom_line(
      data     = emm,
      aes(x = dream, y = response,
          colour = direction, linetype = direction, group = direction),
      linewidth = 0.8,
      position = position_dodge(0.4)
    ) +
    # emmeans means as points on top of lines
    geom_point(
      data     = emm,
      aes(x = dream, y = response, colour = direction),
      size     = 3, shape = 16,
      position = position_dodge(0.4)
    ) +
    scale_colour_manual(values = dir_colors, name = NULL) +
    scale_linetype_manual(values = dir_linetypes, name = NULL) +
    scale_x_discrete(labels = dream_x_labels) +
    scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.04)) +
    labs(x = "Dream recall", y = "GC", title = title) +
    theme_panel
}

p2_1 <- plot_a2(A2_REM_Cz_best, "REM | Cz")
p2_2 <- plot_a2(A2_REM_Pz_best, "REM | Pz")
p2_3 <- plot_a2(A2_N2_Cz_best,  "N2  | Cz")
p2_4 <- plot_a2(A2_N2_Pz_best,  "N2  | Pz")

fig2 <- (p2_1 | p2_2 | p2_3 | p2_4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = "Approach 2: GC by direction and dream recall",
    caption = "Points: individual observations (jittered). Diamonds: group means (emmeans).",
    theme   = theme(
      plot.title   = element_text(face = "bold"),
      plot.caption = element_text(size = 7),
      legend.position = "bottom"
    )
  )

print(fig2)
ggsave("fig2_approach2_mimicPanels_scatter.png", fig2, width = 10, height = 4, dpi = 300)


# =============================================================================
# FIGURE 3 — APPROACH 2 EXPLORATORY
# 4 panels (stage × electrode), each faceted by direction
# x = frequency band, colour = dream group
# =============================================================================

plot_a2x <- function(model, title) {
  
  raw <- get_raw_df(model)
  raw$frequency <- factor(raw$frequency, levels = freq_levels)
  raw$direction <- factor(relabel_dir(raw$direction),
                          levels = c("Brain\u2192Heart", "Heart\u2192Brain"))
  
  emm <- get_emm_df(model, ~ dream * direction * frequency)
  emm$frequency <- factor(emm$frequency, levels = freq_levels)
  emm$direction <- factor(relabel_dir(emm$direction),
                          levels = c("Brain\u2192Heart", "Heart\u2192Brain"))
  
  ggplot() +
    geom_jitter(
      data   = raw,
      aes(x = frequency, y = GC_raw, colour = dream),
      width  = 0.12, height = 0, alpha = 0.6, size = 1.2,
      show.legend = FALSE
    ) +
    geom_point(
      data     = emm,
      aes(x = frequency, y = response, colour = dream),
      size     = 3.5, shape = 18,
      position = position_dodge(0.4)
    ) +
    facet_wrap(~ direction, ncol = 2) +
    scale_colour_manual(values = dream_colors, labels = dream_labels,
                        name = NULL) +
    scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.04)) +
    labs(x = "Frequency band", y = "GC", title = title) +
    theme_panel
}

p3_1 <- plot_a2x(A2x_REM_Cz_best, "REM | Cz")
p3_2 <- plot_a2x(A2x_REM_Pz_best, "REM | Pz")
p3_3 <- plot_a2x(A2x_N2_Cz_best,  "N2  | Cz")
p3_4 <- plot_a2x(A2x_N2_Pz_best,  "N2  | Pz")

fig3 <- (p3_1 / p3_2 / p3_3 / p3_4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = "Approach 2 exploratory: GC by frequency, direction, and dream recall",
    caption = "Points: individual observations (jittered). Diamonds: group means (emmeans).",
    theme   = theme(
      plot.title   = element_text(face = "bold"),
      plot.caption = element_text(size = 7),
      legend.position = "bottom"
    )
  )

print(fig3)
ggsave("fig3_approach2_exploratory_scatter.png", fig3,
       width = 10, height = 14, dpi = 300)
# =============================================================================
#                 DHARMa RESIDUAL DIAGNOSTICS — with testQuantiles
# =============================================================================

cat("\n=================================================================\n")
cat("         DHARMa DIAGNOSTICS                                       \n")
cat("=================================================================\n")
cat("\nTests:\n")
cat(" - testUniformity   : are scaled residuals uniformly distributed?\n")
cat("                      (GLMM equivalent of 'are residuals normal?')\n")
cat(" - testDispersion   : is variance correctly specified on average?\n")
cat(" - testQuantiles    : is variance correctly specified across the\n")
cat("                      range of fitted values? (variance-mean shape)\n")
cat(" - testOutliers     : are there observations the model can't account\n")
cat("                      for at all?\n\n")

run_dharma <- function(model, label) {
  cat("\n---", label, "---\n")
  sim <- simulateResiduals(model, n = 1000)
  print(testUniformity(sim))
  print(testDispersion(sim))
  print(testQuantiles(sim))   # CHANGE: added
  print(testOutliers(sim))
  invisible(sim)
}

dharma_A1_REM_Cz_EEGtoHRV <- run_dharma(A1_REM_Cz_EEGtoHRV_best, "A1_REM_Cz_EEGtoHRV")
dharma_A1_REM_Cz_HRVtoEEG <- run_dharma(A1_REM_Cz_HRVtoEEG_best, "A1_REM_Cz_HRVtoEEG")
dharma_A1_REM_Pz_EEGtoHRV <- run_dharma(A1_REM_Pz_EEGtoHRV_best, "A1_REM_Pz_EEGtoHRV")
dharma_A1_REM_Pz_HRVtoEEG <- run_dharma(A1_REM_Pz_HRVtoEEG_best, "A1_REM_Pz_HRVtoEEG")
dharma_A1_N2_Cz_EEGtoHRV  <- run_dharma(A1_N2_Cz_EEGtoHRV_best,  "A1_N2_Cz_EEGtoHRV")
dharma_A1_N2_Cz_HRVtoEEG  <- run_dharma(A1_N2_Cz_HRVtoEEG_best,  "A1_N2_Cz_HRVtoEEG")
dharma_A1_N2_Pz_EEGtoHRV  <- run_dharma(A1_N2_Pz_EEGtoHRV_best,  "A1_N2_Pz_EEGtoHRV")
dharma_A1_N2_Pz_HRVtoEEG  <- run_dharma(A1_N2_Pz_HRVtoEEG_best,  "A1_N2_Pz_HRVtoEEG")

dharma_A2_REM_Cz <- run_dharma(A2_REM_Cz_best, "A2_REM_Cz")
dharma_A2_REM_Pz <- run_dharma(A2_REM_Pz_best, "A2_REM_Pz")
dharma_A2_N2_Cz  <- run_dharma(A2_N2_Cz_best,  "A2_N2_Cz")
dharma_A2_N2_Pz  <- run_dharma(A2_N2_Pz_best,  "A2_N2_Pz")

dharma_A2x_REM_Cz <- run_dharma(A2x_REM_Cz_best, "A2x_REM_Cz")
dharma_A2x_REM_Pz <- run_dharma(A2x_REM_Pz_best, "A2x_REM_Pz")
dharma_A2x_N2_Cz  <- run_dharma(A2x_N2_Cz_best,  "A2x_N2_Cz")
dharma_A2x_N2_Pz  <- run_dharma(A2x_N2_Pz_best,  "A2x_N2_Pz")

# Save final results
save(pb_results, pb_table, focal_table, nonfocal_table,
     file = "final_results.RData")

cat("\n\nDone. All results saved to final_results.RData\n")

!!!!!!!!!!!!!!!!!!!!!!!!below is old!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# =============================================================================
#                 FOCAL FIXED-EFFECT TESTS (Type III Wald)
# =============================================================================
# Run these first to get the overall pattern. Parametric bootstrap below
# refines the focal terms.

cat("\n=================================================================\n")
cat("         FOCAL FIXED-EFFECT TESTS (Wald chi-squared)             \n")
cat("=================================================================\n")

cat("\n--- Approach 1: A1_REM_Cz_EEGtoHRV ---\n")
print(car::Anova(A1_REM_Cz_EEGtoHRV_best, type = 3))
cat("\n--- Approach 1: A1_REM_Cz_HRVtoEEG ---\n")
print(car::Anova(A1_REM_Cz_HRVtoEEG_best, type = 3))
cat("\n--- Approach 1: A1_REM_Pz_EEGtoHRV ---\n")
print(car::Anova(A1_REM_Pz_EEGtoHRV_best, type = 3))
cat("\n--- Approach 1: A1_REM_Pz_HRVtoEEG ---\n")
print(car::Anova(A1_REM_Pz_HRVtoEEG_best, type = 3))
cat("\n--- Approach 1: A1_N2_Cz_EEGtoHRV ---\n")
print(car::Anova(A1_N2_Cz_EEGtoHRV_best, type = 3))
cat("\n--- Approach 1: A1_N2_Cz_HRVtoEEG ---\n")
print(car::Anova(A1_N2_Cz_HRVtoEEG_best, type = 3))
cat("\n--- Approach 1: A1_N2_Pz_EEGtoHRV ---\n")
print(car::Anova(A1_N2_Pz_EEGtoHRV_best, type = 3))
cat("\n--- Approach 1: A1_N2_Pz_HRVtoEEG ---\n")
print(car::Anova(A1_N2_Pz_HRVtoEEG_best, type = 3))

cat("\n--- Approach 2: A2_REM_Cz ---\n")
print(car::Anova(A2_REM_Cz_best, type = 3))
cat("\n--- Approach 2: A2_REM_Pz ---\n")
print(car::Anova(A2_REM_Pz_best, type = 3))
cat("\n--- Approach 2: A2_N2_Cz ---\n")
print(car::Anova(A2_N2_Cz_best, type = 3))
cat("\n--- Approach 2: A2_N2_Pz ---\n")
print(car::Anova(A2_N2_Pz_best, type = 3))

cat("\n--- Approach 2 exploratory: A2x_REM_Cz ---\n")
print(car::Anova(A2x_REM_Cz_best, type = 3))
cat("\n--- Approach 2 exploratory: A2x_REM_Pz ---\n")
print(car::Anova(A2x_REM_Pz_best, type = 3))
cat("\n--- Approach 2 exploratory: A2x_N2_Cz ---\n")
print(car::Anova(A2x_N2_Cz_best, type = 3))
cat("\n--- Approach 2 exploratory: A2x_N2_Pz ---\n")
print(car::Anova(A2x_N2_Pz_best, type = 3))


# =============================================================================
#                 PARAMETRIC BOOTSTRAP LR TESTS (small-sample valid) => ADJUST SO ITS ON ALL TESTS, NOT JUST FOCAL TESTS!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# =============================================================================
# Tests focal interaction terms via simulation rather than asymptotic chi-sq.
# Slow (~1-2 min per test). Approved by prof as preferred over Wald given the
# sample size and the categorical predictor with >2 levels (frequency).
# Toggle RUN_PBOOT at top of script.

if (RUN_PBOOT) {
  
  cat("\n=================================================================\n")
  cat("         PARAMETRIC BOOTSTRAP LR TESTS                           \n")
  cat("=================================================================\n")
  
  # Approach 1: test dream:frequency interaction (focal term)
  cat("\n--- A1_REM_Cz_EEGtoHRV: dream x frequency ---\n")
  m_full <- A1_REM_Cz_EEGtoHRV_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A1_REM_Cz_HRVtoEEG: dream x frequency ---\n")
  m_full <- A1_REM_Cz_HRVtoEEG_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A1_REM_Pz_EEGtoHRV: dream x frequency ---\n")
  m_full <- A1_REM_Pz_EEGtoHRV_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A1_REM_Pz_HRVtoEEG: dream x frequency ---\n")
  m_full <- A1_REM_Pz_HRVtoEEG_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A1_N2_Cz_EEGtoHRV: dream x frequency ---\n")
  m_full <- A1_N2_Cz_EEGtoHRV_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A1_N2_Cz_HRVtoEEG: dream x frequency ---\n")
  m_full <- A1_N2_Cz_HRVtoEEG_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A1_N2_Pz_EEGtoHRV: dream x frequency ---\n")
  m_full <- A1_N2_Pz_EEGtoHRV_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A1_N2_Pz_HRVtoEEG: dream x frequency ---\n")
  m_full <- A1_N2_Pz_HRVtoEEG_best
  m_red  <- update(m_full, . ~ . - dream:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  # Approach 2: test dream:direction interaction (focal term for hypothesis 3)
  cat("\n--- A2_REM_Cz: dream x direction ---\n")
  m_full <- A2_REM_Cz_best
  m_red  <- update(m_full, . ~ . - dream:direction)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A2_REM_Pz: dream x direction ---\n")
  m_full <- A2_REM_Pz_best
  m_red  <- update(m_full, . ~ . - dream:direction)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A2_N2_Cz: dream x direction ---\n")
  m_full <- A2_N2_Cz_best
  m_red  <- update(m_full, . ~ . - dream:direction)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A2_N2_Pz: dream x direction ---\n")
  m_full <- A2_N2_Pz_best
  m_red  <- update(m_full, . ~ . - dream:direction)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  # Approach 2 exploratory: test the three-way interaction
  cat("\n--- A2x_REM_Cz: dream x direction x frequency ---\n")
  m_full <- A2x_REM_Cz_best
  m_red  <- update(m_full, . ~ . - dream:direction:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A2x_REM_Pz: dream x direction x frequency ---\n")
  m_full <- A2x_REM_Pz_best
  m_red  <- update(m_full, . ~ . - dream:direction:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A2x_N2_Cz: dream x direction x frequency ---\n")
  m_full <- A2x_N2_Cz_best
  m_red  <- update(m_full, . ~ . - dream:direction:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
  
  cat("\n--- A2x_N2_Pz: dream x direction x frequency ---\n")
  m_full <- A2x_N2_Pz_best
  m_red  <- update(m_full, . ~ . - dream:direction:frequency)
  print(PBmodcomp(m_full, m_red, nsim = N_PBOOT))
}


# =============================================================================
#                 FDR CORRECTION on focal p-values
# =============================================================================
# Apply Benjamini-Hochberg jointly across the focal interaction p-values
# from Approach 1 (8 tests of dream x frequency) per prereg approach.

cat("\n=================================================================\n")
cat("         FDR CORRECTION (Benjamini-Hochberg)                     \n")
cat("=================================================================\n")

# Pull dream:frequency p-values from Approach 1 Wald tests
focal_p_A1 <- c(
  A1_REM_Cz_EEGtoHRV = car::Anova(A1_REM_Cz_EEGtoHRV_best, type=3)["dream:frequency", "Pr(>Chisq)"],
  A1_REM_Cz_HRVtoEEG = car::Anova(A1_REM_Cz_HRVtoEEG_best, type=3)["dream:frequency", "Pr(>Chisq)"],
  A1_REM_Pz_EEGtoHRV = car::Anova(A1_REM_Pz_EEGtoHRV_best, type=3)["dream:frequency", "Pr(>Chisq)"],
  A1_REM_Pz_HRVtoEEG = car::Anova(A1_REM_Pz_HRVtoEEG_best, type=3)["dream:frequency", "Pr(>Chisq)"],
  A1_N2_Cz_EEGtoHRV  = car::Anova(A1_N2_Cz_EEGtoHRV_best,  type=3)["dream:frequency", "Pr(>Chisq)"],
  A1_N2_Cz_HRVtoEEG  = car::Anova(A1_N2_Cz_HRVtoEEG_best,  type=3)["dream:frequency", "Pr(>Chisq)"],
  A1_N2_Pz_EEGtoHRV  = car::Anova(A1_N2_Pz_EEGtoHRV_best,  type=3)["dream:frequency", "Pr(>Chisq)"],
  A1_N2_Pz_HRVtoEEG  = car::Anova(A1_N2_Pz_HRVtoEEG_best,  type=3)["dream:frequency", "Pr(>Chisq)"]
)

cat("\nFocal p-values (dream:frequency interaction) from Approach 1:\n")
print(data.frame(
  contrast = names(focal_p_A1),
  p_raw    = focal_p_A1,
  p_FDR    = p.adjust(focal_p_A1, method = "BH")
), row.names = FALSE)


# =============================================================================
#                 EFFECT SIZES: rate ratios via emmeans
# =============================================================================
# exp(beta) gives the multiplicative change in GC for a one-unit change in
# the predictor. type = "response" back-transforms from the log link.

cat("\n=================================================================\n")
cat("         EFFECT SIZES: rate ratios with 95% CIs                  \n")
cat("=================================================================\n")

# Approach 1: dream contrast (yes vs no) per band
cat("\n--- A1_REM_Cz_EEGtoHRV: dream contrast per band ---\n")
print(pairs(emmeans(A1_REM_Cz_EEGtoHRV_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

cat("\n--- A1_REM_Cz_HRVtoEEG: dream contrast per band ---\n")
print(pairs(emmeans(A1_REM_Cz_HRVtoEEG_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

cat("\n--- A1_REM_Pz_EEGtoHRV: dream contrast per band ---\n")
print(pairs(emmeans(A1_REM_Pz_EEGtoHRV_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

cat("\n--- A1_REM_Pz_HRVtoEEG: dream contrast per band ---\n")
print(pairs(emmeans(A1_REM_Pz_HRVtoEEG_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

cat("\n--- A1_N2_Cz_EEGtoHRV: dream contrast per band ---\n")
print(pairs(emmeans(A1_N2_Cz_EEGtoHRV_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

cat("\n--- A1_N2_Cz_HRVtoEEG: dream contrast per band ---\n")
print(pairs(emmeans(A1_N2_Cz_HRVtoEEG_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

cat("\n--- A1_N2_Pz_EEGtoHRV: dream contrast per band ---\n")
print(pairs(emmeans(A1_N2_Pz_EEGtoHRV_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

cat("\n--- A1_N2_Pz_HRVtoEEG: dream contrast per band ---\n")
print(pairs(emmeans(A1_N2_Pz_HRVtoEEG_best, ~ dream | frequency, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))

# Marginal dream contrast (averaged across frequency, per emmeans use approved by prof)
cat("\n--- Approach 1: marginal dream contrast (averaged over frequency) ---\n")
cat("\nA1_REM_Cz_EEGtoHRV:\n")
print(pairs(emmeans(A1_REM_Cz_EEGtoHRV_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))
cat("\nA1_REM_Cz_HRVtoEEG:\n")
print(pairs(emmeans(A1_REM_Cz_HRVtoEEG_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))
cat("\nA1_REM_Pz_EEGtoHRV:\n")
print(pairs(emmeans(A1_REM_Pz_EEGtoHRV_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))
cat("\nA1_REM_Pz_HRVtoEEG:\n")
print(pairs(emmeans(A1_REM_Pz_HRVtoEEG_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))
cat("\nA1_N2_Cz_EEGtoHRV:\n")
print(pairs(emmeans(A1_N2_Cz_EEGtoHRV_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))
cat("\nA1_N2_Cz_HRVtoEEG:\n")
print(pairs(emmeans(A1_N2_Cz_HRVtoEEG_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))
cat("\nA1_N2_Pz_EEGtoHRV:\n")
print(pairs(emmeans(A1_N2_Pz_EEGtoHRV_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))
cat("\nA1_N2_Pz_HRVtoEEG:\n")
print(pairs(emmeans(A1_N2_Pz_HRVtoEEG_best, ~ dream, type = "response"),
            reverse = TRUE, infer = c(TRUE, TRUE)))


# Approach 2: dream x direction interaction
cat("\n--- A2_REM_Cz: dream x direction ---\n")
emm <- emmeans(A2_REM_Cz_best, ~ dream * direction, type = "response")
print(emm)
cat("\nWithin-group direction contrasts:\n")
print(pairs(emm, by = "dream", type = "response"))
cat("\nWithin-direction dream contrasts:\n")
print(pairs(emm, by = "direction", reverse = TRUE, type = "response"))

cat("\n--- A2_REM_Pz: dream x direction ---\n")
emm <- emmeans(A2_REM_Pz_best, ~ dream * direction, type = "response")
print(emm)
cat("\nWithin-group direction contrasts:\n")
print(pairs(emm, by = "dream", type = "response"))
cat("\nWithin-direction dream contrasts:\n")
print(pairs(emm, by = "direction", reverse = TRUE, type = "response"))

cat("\n--- A2_N2_Cz: dream x direction ---\n")
emm <- emmeans(A2_N2_Cz_best, ~ dream * direction, type = "response")
print(emm)
cat("\nWithin-group direction contrasts:\n")
print(pairs(emm, by = "dream", type = "response"))
cat("\nWithin-direction dream contrasts:\n")
print(pairs(emm, by = "direction", reverse = TRUE, type = "response"))

cat("\n--- A2_N2_Pz: dream x direction ---\n")
emm <- emmeans(A2_N2_Pz_best, ~ dream * direction, type = "response")
print(emm)
cat("\nWithin-group direction contrasts:\n")
print(pairs(emm, by = "dream", type = "response"))
cat("\nWithin-direction dream contrasts:\n")
print(pairs(emm, by = "direction", reverse = TRUE, type = "response"))


# =============================================================================
#                 NAKAGAWA R^2 (marginal and conditional)
# =============================================================================

cat("\n=================================================================\n")
cat("         NAKAGAWA R^2                                             \n")
cat("=================================================================\n")

cat("\nApproach 1:\n")
cat("A1_REM_Cz_EEGtoHRV:\n"); print(MuMIn::r.squaredGLMM(A1_REM_Cz_EEGtoHRV_best))
cat("A1_REM_Cz_HRVtoEEG:\n"); print(MuMIn::r.squaredGLMM(A1_REM_Cz_HRVtoEEG_best))
cat("A1_REM_Pz_EEGtoHRV:\n"); print(MuMIn::r.squaredGLMM(A1_REM_Pz_EEGtoHRV_best))
cat("A1_REM_Pz_HRVtoEEG:\n"); print(MuMIn::r.squaredGLMM(A1_REM_Pz_HRVtoEEG_best))
cat("A1_N2_Cz_EEGtoHRV:\n");  print(MuMIn::r.squaredGLMM(A1_N2_Cz_EEGtoHRV_best))
cat("A1_N2_Cz_HRVtoEEG:\n");  print(MuMIn::r.squaredGLMM(A1_N2_Cz_HRVtoEEG_best))
cat("A1_N2_Pz_EEGtoHRV:\n");  print(MuMIn::r.squaredGLMM(A1_N2_Pz_EEGtoHRV_best))
cat("A1_N2_Pz_HRVtoEEG:\n");  print(MuMIn::r.squaredGLMM(A1_N2_Pz_HRVtoEEG_best))

cat("\nApproach 2:\n")
cat("A2_REM_Cz:\n"); print(MuMIn::r.squaredGLMM(A2_REM_Cz_best))
cat("A2_REM_Pz:\n"); print(MuMIn::r.squaredGLMM(A2_REM_Pz_best))
cat("A2_N2_Cz:\n");  print(MuMIn::r.squaredGLMM(A2_N2_Cz_best))
cat("A2_N2_Pz:\n");  print(MuMIn::r.squaredGLMM(A2_N2_Pz_best))

cat("\nApproach 2 exploratory:\n")
cat("A2x_REM_Cz:\n"); print(MuMIn::r.squaredGLMM(A2x_REM_Cz_best))
cat("A2x_REM_Pz:\n"); print(MuMIn::r.squaredGLMM(A2x_REM_Pz_best))
cat("A2x_N2_Cz:\n");  print(MuMIn::r.squaredGLMM(A2x_N2_Cz_best))
cat("A2x_N2_Pz:\n");  print(MuMIn::r.squaredGLMM(A2x_N2_Pz_best))


# =============================================================================
#                 DHARMa RESIDUAL DIAGNOSTICS
# =============================================================================
# Simulated residuals for non-Gaussian models. Look for: KS test ns
# (uniform), dispersion test ns (= 1), outlier test ns. Significant
# deviation = potential model misspecification.

cat("\n=================================================================\n")
cat("         DHARMa DIAGNOSTICS                                       \n")
cat("=================================================================\n")
cat("\nTests:\n")
cat(" - testUniformity   : are scaled residuals uniformly distributed?\n")
cat("                      (GLMM equivalent of 'are residuals normal?')\n")
cat(" - testDispersion   : is variance correctly specified on average?\n")
cat(" - testQuantiles    : is variance correctly specified across the\n")
cat("                      range of fitted values? (variance-mean shape)\n")
cat(" - testOutliers     : are there observations the model can't account\n")
cat("                      for at all?\n\n")
 
run_dharma <- function(model, label) {
  cat("\n---", label, "---\n")
  sim <- simulateResiduals(model, n = 1000)
  print(testUniformity(sim))
  print(testDispersion(sim))
  print(testQuantiles(sim))   # CHANGE: added
  print(testOutliers(sim))
  invisible(sim)
}
 
dharma_A1_REM_Cz_EEGtoHRV <- run_dharma(A1_REM_Cz_EEGtoHRV_best, "A1_REM_Cz_EEGtoHRV")
dharma_A1_REM_Cz_HRVtoEEG <- run_dharma(A1_REM_Cz_HRVtoEEG_best, "A1_REM_Cz_HRVtoEEG")
dharma_A1_REM_Pz_EEGtoHRV <- run_dharma(A1_REM_Pz_EEGtoHRV_best, "A1_REM_Pz_EEGtoHRV")
dharma_A1_REM_Pz_HRVtoEEG <- run_dharma(A1_REM_Pz_HRVtoEEG_best, "A1_REM_Pz_HRVtoEEG")
dharma_A1_N2_Cz_EEGtoHRV  <- run_dharma(A1_N2_Cz_EEGtoHRV_best,  "A1_N2_Cz_EEGtoHRV")
dharma_A1_N2_Cz_HRVtoEEG  <- run_dharma(A1_N2_Cz_HRVtoEEG_best,  "A1_N2_Cz_HRVtoEEG")
dharma_A1_N2_Pz_EEGtoHRV  <- run_dharma(A1_N2_Pz_EEGtoHRV_best,  "A1_N2_Pz_EEGtoHRV")
dharma_A1_N2_Pz_HRVtoEEG  <- run_dharma(A1_N2_Pz_HRVtoEEG_best,  "A1_N2_Pz_HRVtoEEG")
 
dharma_A2_REM_Cz <- run_dharma(A2_REM_Cz_best, "A2_REM_Cz")
dharma_A2_REM_Pz <- run_dharma(A2_REM_Pz_best, "A2_REM_Pz")
dharma_A2_N2_Cz  <- run_dharma(A2_N2_Cz_best,  "A2_N2_Cz")
dharma_A2_N2_Pz  <- run_dharma(A2_N2_Pz_best,  "A2_N2_Pz")
 
dharma_A2x_REM_Cz <- run_dharma(A2x_REM_Cz_best, "A2x_REM_Cz")
dharma_A2x_REM_Pz <- run_dharma(A2x_REM_Pz_best, "A2x_REM_Pz")
dharma_A2x_N2_Cz  <- run_dharma(A2x_N2_Cz_best,  "A2x_N2_Cz")
dharma_A2x_N2_Pz  <- run_dharma(A2x_N2_Pz_best,  "A2x_N2_Pz")
 
# Save final results
save(pb_results, pb_table, focal_table, nonfocal_table,
     file = "final_results.RData")
 
cat("\n\nDone. All results saved to final_results.RData\n")