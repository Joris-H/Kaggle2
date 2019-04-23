library(tidyverse)

##gettubg labels
##########################################################################################
act_labels = read_delim("activity_labels.txt", " ", col_names=F, trim_ws = T) 
act_labels = act_labels %>% dplyr::select(X1,X2)
act_labels


labels <- read_delim("RawData/Train/labels_train.txt", " ", col_names = F)
colnames(labels) <- c('trial', 'userid', 'activity', 'start', 'end')
labels

train_labels = labels
##########################################################################################


##function for entropy 
##########################################################################################
entropy  <- function(x, nbreaks = nclass.Sturges(x)) {
  r = range(x)
  x_binned = findInterval(x, seq(r[1], r[2], len= nbreaks))
  h = tabulate(x_binned, nbins = nbreaks) # fast histogram
  p = h/sum(h)
  -sum(p[p>0] * log(p[p>0]))
}
##########################################################################################



##DEFNING TH FUNCTION 
##########################################################################################
extractTimeDomainFeatures <- 
  function(filename, labels = train_labels) {
    # extract user and experimental run ID's from file name
    username = gsub(".+user(\\d+).+", "\\1", filename)
    expname =  gsub(".+exp(\\d+).+", "\\1", filename)
    
    # import the data from the file
    user01 <- read_delim(filename, " ", col_names = F, progress = TRUE, col_types = "ddd")
    
    # select this users activity labels from the `labels` frame for this experimental run
    user_labels <- 
      labels %>% 
      dplyr::filter(userid==as.numeric(username) & trial == as.numeric(expname)) %>% # select rows pertaining to current signals
      mutate(segment = row_number()) %>%                                             # segment identifies different WALKING_UPSTAIRS etc
      gather(start_end, vec, -trial, -userid, -activity, -segment) %>%               # stack start and end on top of each other
      arrange(vec) %>%                                                               # arrange rows in natural order
      mutate(activity = ifelse(start_end == 'end', NA, activity), activity_id = row_number()) # remove activity label `end` time rows
    
    # add activity labels to each sample
    user <- 
      user01 %>% 
      mutate(sample = row_number()-1) %>%
      mutate(activity_id = findInterval(sample, user_labels$vec)) %>%
      left_join(user_labels, by = "activity_id") 
    
    # split in epochs of 128 samples and compute features per epoch
    usertimedom <- 
      user %>%
      mutate(epoch = sample %/% 128) %>% # epoch = 2.56 sec
      group_by(epoch) %>%
      summarise(
        user_id = username, # added to identify user in data frame rows
        exp_id = expname,   # added to identify experimental run in data frame rows
        activity = names(which.max(table(c("-", activity)))),
        sample = sample[1],
        m1 = mean(X1), 
        m2 = mean(X2),
        m3 = mean(X3),
        sd1 = sd(X1), 
        sd2 = sd(X2),
        sd3 = sd(X3),
        avg_sd = mean(sd1 + sd2 + sd3),
        q1_25 = quantile(X1, .25),
        q2_25 = quantile(X2, .25),
        q3_25 = quantile(X3, .25),
        q1_75 = quantile(X1, .75), 
        q2_75 = quantile(X2, .75),
        q3_75 = quantile(X3, .75),
        skew1 = e1071::skewness(X1),
        skew2 = e1071::skewness(X2),
        skew3 = e1071::skewness(X3),
        AR1.2 = cor(X1, lag(X1, n = 2), use = "pairwise"),
        AR2.2 = cor(X2, lag(X2, n = 2), use = "pairwise"),
        AR3.2 = cor(X3, lag(X3, n = 2), use = "pairwise"),
        AR12.1 = cor(X1, lag(X2), use = "pairwise"),
        AR23.1 = cor(X2, lag(X3), use = "pairwise"),
        AR13.1 = cor(X1, lag(X3), use = "pairwise"),
        max.1 = max(abs(X1)),
        bim.1 = bimodality_coefficient(X1),
        bim.2 = bimodality_coefficient(X2),
        bim.3 = bimodality_coefficient(X3),
        entropy.1 = entropy(X1), 
        entropy.2 = entropy(X2), 
        entropy.3 = entropy(X3), 
        # ... your own features ... (to get inspired, look at the histograms above)
        n=n()
      ) 
    
    usertimedom 
  }
##########################################################################################



###running the function on accelartion and gyro and combining them
##########################################################################################
filenames_acc <- dir("RawData/Train/", "^acc", full.names = TRUE) # for demo only first 5 files

myData_Acc =  filenames_acc %>%
  map_dfr(extractTimeDomainFeatures) # map_dfr runs `extractTimeDomainFeatures` on all elements in filenames and binds results row wise


filenames_gyro <- dir("RawData/Train/", "^gyr", full.names = TRUE) 

myData_Gyro =  filenames_gyro %>%
  map_dfr(extractTimeDomainFeatures) # map_dfr runs `extractTimeDomainFeatures` on all elements in filenames and binds results row wise


myData_Full <- left_join(myData_Acc, myData_Gyro, by = c("epoch", "user_id", "exp_id")) %>% 
  dplyr::select(-activity.y, -sample.y, -n.y) 

myData_Full$activity.x <-  as.integer(myData_Full$activity.x)

##adding the mean differences to the data
myData_Full <- myData_Full %>% mutate(mean_diff_1 = m1.x - lag(m1.x), mean_diff_2 = m2.x - lag(m2.x), mean_diff_3 = m3.x - lag(m3.x)) %>% filter(epoch != 0 )


myData_Full$activity.x <- plyr::mapvalues(myData_Full$activity.x, 1:12, act_labels$X2)
data_to_work_with <- myData_Full  %>% filter(n.x >100) %>% dplyr::select( -c(1,3, sample.x, n.x)) %>% drop_na() %>% mutate(activity.x = factor(activity.x))



##########################################################################################



##This are some remaining issues from Rauls kernel
##########################################################################################
#You should maybe think about the following issues

#What do you do with the unlabelled epochs? That is, the epochs that are marked with -???

#For the competition submision you should only provide predictions for the epochs defined 
#in the sample submission file???

#Should you remove the epochs that do not consist of 128 samples (i.e., n ≠ 128 in the above 
#data frame)????

##########################################################################################

###can use this to split data into training and test based on the users
##########################################################################################
ids <- group_indices(data_to_work_with, user_id)

data_to_work_with_cv <- cbind(ids, data_to_work_with)

data_to_work_with_cv <- as_tibble(data_to_work_with_cv)

train_index <- sample(x = ids, size = 12, replace = F)
train_set <- data_to_work_with_cv %>% filter(ids %in% train_index)
'%ni%' <- Negate('%in%')
test_set <- data_to_work_with_cv %>% filter(ids %ni% train_index)
##########################################################################################



##Function for the Test set
##########################################################################################
filenames_acc_test <- dir("RawData/Test/", "^acc", full.names = TRUE) 
filenames_gyro_test <- dir("RawData/Test/", "^gyr", full.names = TRUE) 



test_features <- function(filenames_acc_test) {
  
  # extract user and experimental run ID's from file name
  username_test = gsub(".+user(\\d+).+", "\\1", filenames_acc_test)
  expname_test =  gsub(".+exp(\\d+).+", "\\1", filenames_acc_test) 
  # import the data from the file
  user01_test <- read_delim(filenames_acc_test, " ", col_names = F)
  
  
  
  user <-  user01_test %>% 
    mutate(sample = row_number()-1) 
  
  usertimedom <- 
    user %>%
    mutate(epoch = sample %/% 128) %>% # epoch = 2.56 sec
    group_by(epoch) %>%
    summarise(
      user_id = username_test, # added to identify user in data frame rows
      exp_id = expname_test,   # added to identify experimental run in data frame rows
      sample = sample[1],
      m1 = mean(X1), 
      m2 = mean(X2),
      m3 = mean(X3),
      sd1 = sd(X1), 
      sd2 = sd(X2),
      sd3 = sd(X3),
      avg_sd = mean(sd1 + sd2 + sd3),
      q1_25 = quantile(X1, .25),
      q2_25 = quantile(X2, .25),
      q3_25 = quantile(X3, .25),
      q1_75 = quantile(X1, .75), 
      q2_75 = quantile(X2, .75),
      q3_75 = quantile(X3, .75),
      skew1 = e1071::skewness(X1),
      skew2 = e1071::skewness(X2),
      skew3 = e1071::skewness(X3),
      AR1.2 = cor(X1, lag(X1, n = 2), use = "pairwise"),
      AR2.2 = cor(X2, lag(X2, n = 2), use = "pairwise"),
      AR3.2 = cor(X3, lag(X3, n = 2), use = "pairwise"),
      AR12.1 = cor(X1, lag(X2), use = "pairwise"),
      AR23.1 = cor(X2, lag(X3), use = "pairwise"),
      AR13.1 = cor(X1, lag(X3), use = "pairwise"),
      sum.1 = sum(X1),
      sum.2 = sum(X2),
      sum.3 = sum(X3),
      max.1 = max(abs(X1)),
      bim.1 = bimodality_coefficient(X1),
      bim.2 = bimodality_coefficient(X2),
      bim.3 = bimodality_coefficient(X3),
      entropy.1 = entropy(X1), 
      entropy.2 = entropy(X2), 
      entropy.3 = entropy(X3), 
      # ... your own features ... (to get inspired, look at the histograms above)
      n=n()
    ) 
  usertimedom
}
##########################################################################################


##Making the test datast
##########################################################################################
myData_Acc_test =  filenames_acc_test %>%
  map_dfr(test_features)

myData_Gryo_test <- filenames_gyro_test %>%
  map_dfr(test_features)

myData_Full_test <- left_join(myData_Acc_test, myData_Gryo_test, by = c("epoch", "user_id", "exp_id")) %>% 
  dplyr::select(-sample.y, -n.y, -n.x) 

##adding mean differences to the test data set
myData_Full_test <- myData_Full_test %>% mutate(mean_diff_1 = m1.x - lag(m1.x), mean_diff_2 = m2.x - lag(m2.x), mean_diff_3 = m3.x - lag(m3.x)) %>% filter(epoch != 0 )

myData_Full_test %>%
  mutate(user_id = paste("user", user_id, sep=""), exp_id = paste("exp", exp_id, sep="")) %>%
  unite(Id, user_id, exp_id, sample.x) %>%
  dplyr::select(Id, Predicted = activity) %>%
  write_csv("test_set_predictions.csv")

file.show("test_set_predictions.csv")


