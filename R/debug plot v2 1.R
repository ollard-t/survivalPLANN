# path1 <- "~/Documents/Rstudio/Simulations/BASES vFinale/"
path1 <- "~/Documents/Rstudio/Simulations/BASES/"

ind = sample(1:1000, 8, replace = FALSE)
for(i in ind){
data <- read.csv(file = paste0(path1,"ind3000NPH_HC/",i,"_df3000.csv"), sep = ";")
data_train <- data
newtimes = seq(1, 3652.41, by = 365.241/12)
data_train_HC <- data_train[data_train$sex == 1 & data_train$colon == 1,]
data_train_HR <- data_train[data_train$sex == 1 & data_train$colon == 0,]

data_train_FC <- data_train[data_train$sex == 2 & data_train$colon == 1,]
data_train_FR <- data_train[data_train$sex == 2 & data_train$colon == 0,]

strata_names = c("HC", "HR", "FC", "FR")


flexmodel <- survivalFLEXNET2(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + 
      strata(sex.organ) + ratetable(age, year, sexchara), data = data_train,
      ratetable=slopop, m = 2, mpos = NULL, mquant = NULL, init = NULL, delta_th = 0, weights=NULL)

pred <- predict(flexmodel, newtimes = newtimes)
# ind = as.numeric(row.names(data_train_HC))
# HC.mean <- apply(predictions[ind,], mean, MARGIN = 2)
# 
# ind = as.numeric(row.names(data_train_HR))
# HR.mean <- apply(predictions[ind,], mean, MARGIN = 2)
# 
# ind = as.numeric(row.names(data_train_FC))
# FC.mean <- apply(predictions[ind,], mean, MARGIN = 2)
# 
# ind = as.numeric(row.names(data_train_FR))
# FR.mean <- apply(predictions[ind,], mean, MARGIN = 2)
#####################
for(j in strata_names){
  
  data_train_var <- get(paste0("data_train_", j))
  
  PPS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train_var,
                          ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                          rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  poharpredS <- PPS$surv
  
  assign(paste0("PP_", j), PPS)
  assign(paste0("poharpred_", j), poharpredS)
}

for(j in strata_names){
  
  data_train_var <- get(paste0("data_train_", j))
  
  flexpred1.2S <- predict(flexmodel, newtimes = newtimes, newdata = data_train_var)$predictions
  
  mean.flex1.2S <- apply(flexpred1.2S, FUN="mean", MARGIN=2)
  
  assign(paste0("flexpred1.2_", j), flexpred1.2S)
  assign(paste0("mean.flex1.2_", j), mean.flex1.2S)
}

#####################
plot(newtimes, poharpred_HC, type = 'l', xlim = c(0,3650), ylim = c(0,1), main = i)
lines(newtimes, poharpred_HR, col = 'red')
lines(newtimes, poharpred_FC, col = 'green')
lines(newtimes, poharpred_FR, col = 'blue')

lines(newtimes, mean.flex1.2_HC, col = 'black')
lines(newtimes, mean.flex1.2_HR, col = 'red')
lines(newtimes, mean.flex1.2_FC, col = 'green')
lines(newtimes, mean.flex1.2_FR, col = 'blue')
}
