rm(list=ls()) #clear workspace

#load packages
library(tidyverse)      #for tidying and wrangling
library(dplyr)          #useful tidying and wrangling functions

#load data
ANON <- read_csv("exp4_anon.csv")

### ---- ESTIMATING EFFECT OF LENGTH AND POSITION ON RT ---- ###

FILLERS <- ANON %>%
  #filter to leave only FILL trials
  filter(StimType == "fill") %>%
  #long format with POS(ITION) column from CHUNK1, CHUNK2 etc.
  pivot_longer(cols = c(starts_with("CHUNK"), "SPILL"),
               names_to = "POS",
               values_to = "CHUNK",
               names_pattern = "CHUNK(\\d+)|SPILL",
               names_transform = list(POS = as.integer)) %>%
  mutate(POS = case_when(is.na(POS) ~ "SPILL",
                         TRUE ~ as.character(POS))) %>%
  #group by Participant as want Participant's RT
  #for given chunk
  group_by(Participant) %>%
  #populate RT column with given chunk RT
  mutate(RT = case_when(POS == 1 ~ C1.RT,
                        POS == 2 ~ C2.RT,
                        POS == 3 ~ C3.RT,
                        POS == 4 ~ C4.RT,
                        POS == "SPILL" ~ SPILL.RT)) %>%
  #add LEN(GTH) column for chunk length
  mutate(LEN = str_count(CHUNK, '\\S+'))

EXP <- ANON %>%
  #filter to leave only EXP trials
  filter(StimType == "exp") %>%
  #long format with POS(ITION) column from CHUNK1, CHUNK2 etc.
  pivot_longer(cols = c(starts_with("CHUNK"), "SPILL"),
               names_to = "POS",
               values_to = "CHUNK",
               names_pattern = "CHUNK(\\d+)|SPILL",
               names_transform = list(POS = as.integer)) %>%
  mutate(POS = case_when(is.na(POS) ~ "SPILL",
                         TRUE ~ as.character(POS))) %>%
  #group by Participant as want Participant's RT
  #for given chunk
  group_by(Participant) %>%
  #populate RT column with given chunk RT
  mutate(RT = case_when(POS == 1 ~ C1.RT,
                        POS == 2 ~ C2.RT,
                        POS == 3 ~ C3.RT,
                        POS == 4 ~ C4.RT,
                        POS == "SPILL" ~ SPILL.RT)) %>%
  #add LEN(GTH) column for chunk length
  mutate(LEN = str_count(CHUNK, '\\S+'))

### ---- GENERATE SOME COUNTS ---- ###

#get all LEN + POS combos
all_combinations <- expand.grid(
  LEN = unique(FILLERS$LEN),
  POS = unique(FILLERS$POS),
  stringsAsFactors = FALSE
)

#get actual LEN + POS counts
length_position_counts <- FILLERS %>%
  group_by(LEN, POS) %>%
  summarise(Count = n_distinct(CHUNK), .groups = "drop")

exc_table <- all_combinations %>%
  filter(LEN %in% c(4, 5),
         POS %in% c(3, 4)) %>%
  left_join(length_position_counts, by = c("LEN", "POS")) %>%
  mutate(Count = replace_na(Count, 0)) %>%
  arrange(LEN) %>%
  xtable()

all_table <- all_combinations %>%
  left_join(length_position_counts, by = c("LEN", "POS")) %>%
  mutate(Count = replace_na(Count, 0)) %>%
  pivot_wider(names_from = POS, values_from = Count, values_fill = 0) %>%
  arrange(LEN) %>%
  xtable()

### ---- CALC PREDICTED RTS FOR EACH PARTICIPANT BASED ON LEN AND POS **INDEPENDENTLY** ---- ###

#new data frame
final_predictions <- EXP

#new columns in to store predicted RTs
final_predictions$LEN.PRED <- NA
final_predictions$POS.PRED <- NA

#loop through each participant
for(participant in unique(FILLERS$Participant)) {
  # Filter data for the current participant in the filler dataset
  participant_data <- FILLERS %>% filter(Participant == participant)
  
  #fit individual model
  len_mod <- lm(RT ~ LEN, data = participant_data)
  pos_mod <- lm(RT ~ POS, data = participant_data)

  #filter prediction data for current participant in long dataset
  prediction_data <- EXP %>% filter(Participant == participant)
  
  #predict RT
  len_preds <- predict(len_mod, newdata = prediction_data)
  pos_preds <- predict(pos_mod, newdata = prediction_data)

  #append predictions in final_predictions data frame
  final_predictions$LEN.PRED[final_predictions$Participant == participant] <- len_preds
  final_predictions$POS.PRED[final_predictions$Participant == participant] <- pos_preds

}

#now need to convert back to wide format to export for analysis
#wrangle so we have a predicted RT for each trial for each participant
#for each of the four chunks

PREDS_LEN <- final_predictions %>%
  dplyr::select(Participant, TrialN, POS, LEN.PRED) %>%
  pivot_wider(names_from = POS,
              values_from = LEN.PRED,
              names_glue = "C{POS}.LEN.PRED") %>%
  rename(SPILL.LEN.PRED = CSPILL.LEN.PRED)

PREDS_POS <- final_predictions %>%
  dplyr::select(Participant, TrialN, POS, POS.PRED) %>%
  pivot_wider(names_from = POS,
              values_from = POS.PRED,
              names_glue = "C{POS}.POS.PRED") %>%
  rename(SPILL.POS.PRED = CSPILL.POS.PRED)

#join the predicted RTs and the experimental data
exp4_tidied <- ANON %>%
  filter(StimType == "exp") %>%
  left_join(PREDS_LEN, by = c("Participant", "TrialN")) %>%
  left_join(PREDS_POS, by = c("Participant", "TrialN"))

### ---- VP RESIDUALS SINCE THERE ARE ENOUGH OBS FOR THAT ---- ###

final_predictions$VP.PRED <- NA   # Create a column for VP predictions

for (participant in unique(FILLERS$Participant)) {
  
  participant_data <- FILLERS %>% filter(Participant == participant)
  
  vp_mod <- lm(RT ~ LEN + POS, data = participant_data)   # Model for VP.RT
  
  vp_data <- EXP %>% filter(Participant == participant)
  
  vp_preds <- predict(vp_mod, newdata = vp_data)
  
  final_predictions$VP.PRED[final_predictions$Participant == participant] <- vp_preds
}

#wide format
PREDS_VP <- final_predictions %>%
  dplyr::select(Participant, TrialN, POS, VP.PRED) %>%
  pivot_wider(names_from = POS,
              values_from = VP.PRED,
              names_glue = "C{POS}.VP.PRED") %>%
  dplyr::select(Participant, TrialN, C3.VP.PRED, C4.VP.PRED)

#join predictions
exp4_tidied <- exp4_tidied %>%
  filter(StimType == "exp") %>%
  left_join(PREDS_VP, by = c("Participant", "TrialN"))

### ---- WRANGLING AND TIDYING ---- ###

#working out cumulative counts for info source - repetition effect
CumulativeCounts <- exp4_tidied %>%
  group_by(Participant, MAXRestrType) %>%  #group by participant and restr type
  arrange(Participant, TrialN) %>%  #ensure order is correct
  mutate(MAXRestrTypeCount = cumsum(MAXRestrType == MAXRestrType)) %>%  #cumulative count calc.
  pivot_wider(names_from = MAXRestr, values_from = MAXRestrTypeCount, 
              id_cols = c(Participant, TrialN, MAXRestrType)) %>% 
  mutate(InfoCount = case_when(
    MAXRestrType == "vlog" ~ `in the vlog`,
    MAXRestrType == "brochure" ~ `in the brochure`,
    MAXRestrType == "app" ~ `on the app`,
    MAXRestrType == "blog" ~ `in the blog`,
    MAXRestrType == "website" ~ `on the website`
  )) %>%
  dplyr::select(Participant, TrialN, MAXRestrType, InfoCount)

#adding cumulative counts
exp4_tidied <- exp4_tidied %>%
  right_join(CumulativeCounts, by = c("Participant", "TrialN", "MAXRestrType")) %>%
  dplyr::select(-c(MENRestr, MENRestrType, MAXRestr,
                   MAXRestrType, StimType)) #remove irrelevant columns

#general tidying/wrangling
exp4_tidied <- exp4_tidied %>%
  mutate(ExceptiveType = case_when(ExceptiveType == "MEN" ~ "SUB",
                                    ExceptiveType == "MAX" ~ "MAX")) %>%
  mutate(Condition = str_c(ExceptivePosition, "-", ExceptiveType)) %>%
  mutate(across(ends_with(".RT"), log, .names = "{.col}.LOG")) %>% #logged RTs
  mutate(Condition = str_c(ExceptivePosition, "-", ExceptiveType)) %>%
  mutate(EXC.RT = case_when(ExceptivePosition == "PRE" ~ C3.RT,
                           ExceptivePosition == "PST" ~ C4.RT),
         EXC.LOG = case_when(ExceptivePosition == "PRE" ~ C3.RT.LOG,
                            ExceptivePosition == "PST" ~ C4.RT.LOG),
         VP.RT = case_when(ExceptivePosition == "PRE" ~ C4.RT,
                           ExceptivePosition == "PST" ~ C3.RT),
         VP.LOG = case_when(ExceptivePosition == "PRE" ~ C4.RT.LOG,
                            ExceptivePosition == "PST" ~ C3.RT.LOG)) %>%
  #predicted reading time columns
  mutate(EXC.LEN = case_when(ExceptivePosition == "PRE" ~ C3.LEN.PRED,
                             ExceptivePosition == "PST" ~ C4.LEN.PRED),
         EXC.POS = case_when(ExceptivePosition == "PRE" ~ C3.POS.PRED,
                             ExceptivePosition == "PST" ~ C4.POS.PRED)) %>%
  mutate(VP.PRED = case_when(ExceptivePosition == "PST" ~ C3.VP.PRED,
                      ExceptivePosition == "PRE" ~ C4.VP.PRED)) %>%
  mutate(VP.RESID = VP.RT - VP.PRED) #calc VP.RESID value

#export dataset
exp4_tidied %>%
  write_csv("exp4_tidied.csv")

#mean age of participants for write-up
exp4_tidied %>%
  summarise(MeanAge = mean(Age, na.rm = TRUE),
            SDAge = sd(Age, na.rm = TRUE))