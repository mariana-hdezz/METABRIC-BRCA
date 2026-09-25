library(rsample)

metadata_ER_POS_REC <- readRDS("./output_data/metadata_ER_POS_REC.RDS")
metadata_ER_POS_SURV <- readRDS("./output_data/metadata_ER_POS_SURV.RDS")

# Survival split ----------------------------------------------------------


dim(metadata_ER_POS_SURV)


set.seed(123)

metabric_surv_split <- initial_split(
  metadata_ER_POS_SURV,
  prop = 0.8,
  strata = EVENT_STAT)

metadata_surv_train <- training(metabric_surv_split)

metadata_surv_test  <- testing(metabric_surv_split)

train_surv_id <- metadata_surv_train$PATIENT_ID
test_surv_id <- metadata_surv_test$PATIENT_ID

length(train_surv_id)
length(test_surv_id)

table(metadata_surv_train$EVENT_STAT)
table(metadata_surv_test$EVENT_STAT) 

intersect(train_surv_id, test_surv_id)


# Recurrence --------------------------------------------------------------


dim(metadata_ER_POS_REC)


set.seed(123)

metabric_rec_split <- initial_split(
  metadata_ER_POS_REC,
  prop = 0.8,
  strata = EVENT_STAT)

metadata_rec_train <- training(metabric_rec_split)

metadata_rec_test  <- testing(metabric_rec_split)

train_rec_id <- metadata_rec_train$PATIENT_ID
test_rec_id <- metadata_rec_test$PATIENT_ID

length(train_rec_id)
length(test_rec_id)



saveRDS(train_rec_id, "./output_data/train_rec_id.RDS")
saveRDS(test_rec_id, "./output_data/test_rec_id.RDS")
saveRDS(metadata_rec_train, "./output_data/metadata_rec_train.RDS")
saveRDS(metadata_rec_test,  "./output_data/metadata_rec_test.RDS")

saveRDS(train_surv_id, "./output_data/train_surv_id.RDS")
saveRDS(test_surv_id, "./output_data/test_surv_id.RDS")
saveRDS(metadata_surv_train, "./output_data/metadata_surv_train.RDS")
saveRDS(metadata_surv_test,  "./output_data/metadata_surv_test.RDS")


rm(list = ls())
gc()
