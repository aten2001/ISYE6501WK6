library('glmnet')
library('glmnetUtils')

q2_data <- read.csv('breast-cancer-wisconsin.data', header = FALSE, na.strings = '?')

# for (i in 1:ncol(q2_data)) {
#     q2_data[,i] <- factor(q2_data[,i])
# }

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

impute_mode <- function(data) {
    # Impute using mode (col 7 is categorical)
    ind.missing <- is.na(data[,7])

    mode <- Mode(data[,7])

    data.imp <- data
    data.imp[ind.missing,7] <- mode
    
    return(data.imp)

}

impute_reg <- function(data) {
    # Impute using LASSO regression
    ind <- is.na(data[,7])
    
    data.tbi <- data[ind, ]
    data.train <- data[!ind, ]
    
    model <- cv.glmnet(V7 ~., data = data.train, family = 'gaussian')
    
    pred <- predict(model, newdata = data.tbi, type = 'response', s = 'lambda.min')
    
    pred <- floor(0.5 + pred)
    
    data.imp <- data
    data.imp[ind,7] <- pred
    data.imp[,7] <- factor(data.imp[,7])
    
    return(data.imp)
}

impute_perturb <- function(data) {
    # Impute using LASSO regression
    ind <- is.na(data[,7])
    
    data.tbi <- data[ind, ]
    data.train <- data[!ind, ]
    
    model <- cv.glmnet(V7 ~., data = data.train, family = 'gaussian')
    
    pred <- predict(model, newdata = data.tbi, type = 'response', s = 'lambda.min')
    
    noise <- rnorm(n = length(pred), mean = 0, 
                   sd = sd(as.numeric(as.character(data.train[,7]))))
    
    pred <- min(max(1, floor(pred + noise + 0.5)), 10)
    
    data.imp <- data
    data.imp[ind,7] <- pred
    data.imp[,7] <- factor(data.imp[,7])
    
    return(data.imp)
}

set.seed(123)
q2.train.imp.mode <- impute_mode(q2_data.train)
q2.train.imp.reg <- impute_reg(q2_data.train)
q2.train.imp.pert <- impute_perturb(q2_data.train)