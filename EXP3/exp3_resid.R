rm(list=ls()) #clear workspace

#load packages
library(tidyverse)      #for tidying and wrangling
library(dplyr)          #useful tidying and wrangling functions

#load data
ANON <- read_csv("exp3_anon.csv")

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
  mutate(LEN = str_count(CHUNK, '\\w+')) %>%
  mutate(LOG = log(RT))

#do a long format version of exp6_tidied_EXP same as above
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
  mutate(LEN = str_count(CHUNK, '\\w+')) %>%
  mutate(LOG = log(RT))

#now to calculate the predicted RTs for each participant,
#for RAW, LOG,

#new data frame
final_predictions <- EXP

#new columns in to store predicted RTs
final_predictions$PRED.RAW <- NA
final_predictions$PRED.LOG <- NA

#loop through each participant
for(participant in unique(FILLERS$Participant)) {
  # Filter data for the current participant in the filler dataset
  participant_data <- FILLERS %>% filter(Participant == participant)
  
  #fit individual model
  raw_mod <- lm(RT ~ LEN + POS, data = participant_data)
  log_mod <- lm(LOG ~ LEN + POS, data = participant_data)

  #filter prediction data for current participant in long dataset
  prediction_data <- EXP %>% filter(Participant == participant)
  
  #predict RT
  raw_preds <- predict(raw_mod, newdata = prediction_data)
  log_preds <- predict(log_mod, newdata = prediction_data)

  #append predictions in final_predictions data frame
  final_predictions$PRED.RAW[final_predictions$Participant == participant] <- raw_preds
  final_predictions$PRED.LOG[final_predictions$Participant == participant] <- log_preds
}

#now need to convert back to wide format to export for analysis
#wrangle so we have a predicted RT for each trial for each participant
#for each of the four chunks

PREDS_RAW <- final_predictions %>%
  dplyr::select(Participant, TrialN, POS, PRED.RAW) %>%
  pivot_wider(names_from = POS,
              values_from = PRED.RAW,
              names_glue = "C{POS}.PRED") %>%
  rename(SPILL.PRED = CSPILL.PRED)

PREDS_LOG <- final_predictions %>%
  dplyr::select(Participant, TrialN, POS, PRED.LOG) %>%
  pivot_wider(names_from = POS,
              values_from = PRED.LOG,
              names_glue = "C{POS}.PRED.LOG") %>%
  rename(SPILL.PRED.LOG = CSPILL.PRED.LOG)

#join the predicted RTs and the experimental data
exp3_tidied <- ANON %>%
  filter(StimType == "exp") %>%
  left_join(PREDS_RAW, by = c("Participant", "TrialN")) %>%
  left_join(PREDS_LOG, by = c("Participant", "TrialN"))

### ---- WRANGLING AND TIDYING ---- ###

exp3_tidied <- exp3_tidied %>%
  mutate(RestrictorType = case_when(RestrictorType == "EM" ~ "SUB",
                                    RestrictorType == "NEM" ~ "OTH")) %>%
  mutate(Condition = str_c(RestrictorPosition, "-", RestrictorType)) %>%
  mutate(across(ends_with(".RT"), log, .names = "{.col}.LOG")) %>% #loged RTs
  mutate(across(ends_with(".RT"), 
                ~ .x - get(sub("\\.RT$", ".PRED", cur_column())), 
                .names = "{.col}.RESID")) %>%
  mutate(across(ends_with(".RT.LOG"), 
                ~ .x - get(sub("\\.RT.LOG$", ".PRED.LOG", cur_column())), 
                .names = "{.col}.RESID")) %>%
  mutate(Condition = str_c(RestrictorPosition, "-", RestrictorType)) %>%
  mutate(PP.RT = case_when(RestrictorPosition == "INT" ~ C3.RT,
                           RestrictorPosition == "EXT" ~ C4.RT),
         PP.LOG = case_when(RestrictorPosition == "INT" ~ C3.RT.LOG,
                            RestrictorPosition == "EXT" ~ C4.RT.LOG),
         VP.RT = case_when(RestrictorPosition == "INT" ~ C4.RT,
                           RestrictorPosition == "EXT" ~ C3.RT),
         VP.LOG = case_when(RestrictorPosition == "INT" ~ C4.RT.LOG,
                            RestrictorPosition == "EXT" ~ C3.RT.LOG)) %>%
  #predicted reading time columns
  mutate(PP.PRED = case_when(RestrictorPosition == "INT" ~ C3.PRED,
                             RestrictorPosition == "EXT" ~ C4.PRED),
         PP.PRED.LOG = case_when(RestrictorPosition == "INT" ~ C3.PRED.LOG,
                                 RestrictorPosition == "EXT" ~ C4.PRED.LOG),
         VP.PRED = case_when(RestrictorPosition == "EXT" ~ C3.PRED,
                             RestrictorPosition == "INT" ~ C4.PRED),
         VP.PRED.LOG = case_when(RestrictorPosition == "EXT" ~ C3.PRED.LOG,
                                 RestrictorPosition == "INT" ~ C4.PRED.LOG)) %>%
  #resid columns
  mutate(PP.RESID = PP.RT - PP.PRED,
         PP.RESID.LOG = PP.LOG - PP.PRED.LOG,
         VP.RESID = VP.RT - VP.PRED,
         VP.RESID.LOG = VP.LOG - VP.PRED.LOG,
         SPILL.RESID = SPILL.RT - SPILL.PRED,
         SPILL.RESID.LOG = SPILL.RT.LOG - SPILL.PRED.LOG)
  
#export dataset
exp3_tidied %>%
  write_csv("exp3_tidied.csv")

#mean age of participants for write-up
exp3_tidied %>%
  summarise(MeanAge = mean(Age, na.rm = TRUE),
            SDAge = sd(Age, na.rm = TRUE))









