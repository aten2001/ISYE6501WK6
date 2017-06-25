library('glmnet')
library('glmnetUtils')
library('caret')
library('parallel')
library('doParallel')

q2_data <- read.csv('breast-cancer-wisconsin.data', header = FALSE, na.strings = '?')

for (i in 1:ncol(q2_data)) {
    q2_data[,i] <- factor(q2_data[,i])
}

# validation/training split (20/80)
set.seed(123)
inTrain <- sample(1:nrow(q2_data), size = ceiling(nrow(q2_data)*0.8))
q2_data.train <- q2_data[inTrain,]
q2_data.validation <- q2_data[-inTrain,]

sum(is.na(q2_data))

percentages <- rep(0, ncol(q2_data))

for (i in 1:ncol(q2_data)) {
    percentages[i] <- sum(is.na(q2_data[,i]))/nrow(q2_data)
}


# Column 7 (Bare Nuclei) contains 16 missing data
# around 2.3% missing


Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}


train_impute_mode <- function(data) {
    mode <- Mode(data[,7])
    
    return(mode)
}

predict_impute_mode <- function(model, newdata) {
    ind.missing <- is.na(newdata[,7])
    
    newdata[ind.missing,7] <- model
    
    return(newdata)
}

train_impute_reg <- function(data) {
    # Impute using LASSO regression
    ind <- is.na(data[,7])

    data.train <- data[!ind, ]
    
    model <- cv.glmnet(V7 ~ V1 + V2 + V3 + V4 + V5 + V6 + V8 + V9 + V10, 
                       data = data.train, family = 'gaussian')
    
    return(model)
}

predict_impute_reg <- function(model, newdata) {
    ind <- is.na(newdata[,7])
    
    data.tbi <- newdata[ind,]
    
    pred <- predict(model, newdata = data.tbi, type = 'response', s = 'lambda.min')
    
    pred <- floor(0.5 + pred)
    
    data.imp <- newdata
    data.imp[ind,7] <- pred
    data.imp[,7] <- factor(data.imp[,7])
    
    return(data.imp)
}


train_impute_perturb <- function(data) {
    # Impute using LASSO regression
    ind <- is.na(data[,7])

    data.train <- data[!ind, ]
    
    model <- cv.glmnet(V7 ~ V1 + V2 + V3 + V4 + V5 + V6 + V8 + V9 + V10, 
                       data = data.train, family = 'gaussian')
    
    return(model)
}

predict_impute_perturb <- function(model, newdata) {
    # Impute with model prediction +  normally distributed noise. 
    # SD = model RMSE.
    ind <- is.na(newdata[,7])
    
    data.tbi <- newdata[ind, ]
    
    sd = sqrt(min(model$cvm))
    
    pred <- predict(model, newdata = data.tbi, type = 'response', s = 'lambda.min')
    
    noise <- rnorm(n = length(pred), mean = 0, 
                   sd = sd)
    
    pred <- pmin(pmax(1, floor(pred + noise + 0.5)), 10)
    
    data.imp <- newdata
    data.imp[ind,7] <- pred
    data.imp[,7] <- factor(data.imp[,7])
    
    return(data.imp)
}


imp.model.mode <- train_impute_mode(q2_data.train)
imp.model.reg <- train_impute_reg(q2_data.train)
imp.model.perturb <- train_impute_perturb(q2_data.train)

q2.train.imp.mode <- predict_impute_mode(imp.model.mode, q2_data.train)
q2.train.imp.reg <- predict_impute_reg(imp.model.reg, q2_data.train)
set.seed(123)
q2.train.imp.pert <- predict_impute_perturb(imp.model.perturb, q2_data.train)

## Final part
missing <- is.na(q2_data.train[,7])
q2.train.complete <- q2_data.train[!missing,]

# Build models & compare.
trCon <- trainControl(method = 'cv')

set.seed(123)

cl <- makeCluster(detectCores(logical = FALSE))
registerDoParallel(cl)
# 
# model.mode <- train(V11 ~., data = q2.train.imp.mode, method = 'rf', trControl = trCon)
# model.reg <- train(V11 ~., data = q2.train.imp.reg, method = 'rf', trControl = trCOn)
# model.pert <- train()