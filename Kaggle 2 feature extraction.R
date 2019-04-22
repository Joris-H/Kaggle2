library(tidyverse)

act_labels = read_delim("activity_labels.txt", " ", col_names=F, trim_ws = T) 
act_labels = act_labels %>% select(X1,X2)
act_labels

labels <- read_delim("RawData/Train/labels_train.txt", " ", col_names = F)
colnames(labels) <- c('trial', 'userid', 'activity', 'start', 'end')
labels


# identify the file name and extract the 'username' (participant ID) and 'expname' (experimental run)
filename = "RawData/Train/acc_exp01_user01.txt"
username = gsub(".+user(\\d+).+", "\\1", filename)
expname =  gsub(".+exp(\\d+).+", "\\1", filename)

# import the data from the file
user01 <- read_delim(filename, " ", col_names = F)

head(user01)

#Each column is a signal. Subsequent rows are subsequent measurement samples, and 
#so we treat rownumber as a time indicator (to keep the distinction clear we'll talk about 
#sample number).
#Let's have a look at the signal wave forms:

plot.ts(user01, xlab="Sample number")

labels[1:2,]


findInterval(.4, c(0,.25,.5,1))

user_labels <-  labels %>% 
  dplyr::filter(userid==as.numeric(username) & trial == as.numeric(expname)) %>% 
  # select rows pertaining to current signals
  mutate(segment = row_number()) %>%                                             
  # segment identifies different WALKING_UPSTAIRS etc (handy)
  gather(start_end, vec, -trial, -userid, -activity, -segment) %>%               
  # stack start and end on top of each other
  arrange(vec) %>%                                                              
  # arrange rows in natural order
  mutate(activity = ifelse(start_end == 'end', NA, activity), activity_id = row_number()) 
# remove activity label from rows with `end` times

user <-  user01 %>% 
  mutate(sample = row_number()-1) %>%
  mutate(activity_id = findInterval(sample, user_labels$vec)) %>%
  left_join(user_labels)

user

user %>% #drop_na() %>%
  ggplot(aes(sample, X1, col = factor(activity), group=segment)) + # first without 'segment', then solve issues
  geom_line()   #+ geom_line(aes(y=X2)) + geom_line(aes(y=X3)) 
#facet_wrap(~segment+activity, scales="free_x")


user %>%
  ggplot(aes(X1)) + 
  geom_histogram(bins=40, fill=1, alpha=0.5) + 
  #geom_histogram(aes(X2), bins=40, fill = 2, alpha=0.5) + geom_histogram(aes(X3), bins=40, fill = 4, alpha=0.5) +
  facet_wrap(~activity, scales = "free_y")


usertimedom <- user %>%
  mutate(epoch = sample %/% 128) %>% # epoch = 2.56 sec
  group_by(epoch) %>%
  summarise(
    activity = names(which.max(table(c("-", activity)))),
    sample = sample[1],
    m1 = mean(X1), 
    m2 = mean(X2), 
    sd1 = sd(X1), 
    q1_25 = quantile(X1, .25),
    skew1 = e1071::skewness(X1),
    AR1.1 = cor(X1, lag(X1), use = "pairwise"),
    AR1.2 = cor(X1, lag(X1, n = 2), use = "pairwise"),
    AR12.1 = cor(X1, lag(X2), use = "pairwise"),
    AutoReg.1 = ar(X1, aic=FALSE, order.max=1)$ar,
    # ... your own features ... (to get inspired, look at the histograms above)
    n=n())

head(usertimedom)
tail(usertimedom)


##DEFNING TH FUNCTION 
train_labels = labels

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
        
        # ... your own features ... (to get inspired, look at the histograms above)
        n=n()
      ) 
    
    usertimedom 
  }

filenames_acc <- dir("RawData/Train/", "^acc", full.names = TRUE) # for demo only first 5 files

myData_Acc =  filenames_acc %>%
  map_dfr(extractTimeDomainFeatures) # map_dfr runs `extractTimeDomainFeatures` on all elements in filenames and binds results row wise

head(myData_Acc)


filenames_gyro <- dir("RawData/Train/", "^gyr", full.names = TRUE) 

myData_Gyro =  filenames_gyro %>%
  map_dfr(extractTimeDomainFeatures) # map_dfr runs `extractTimeDomainFeatures` on all elements in filenames and binds results row wise

head(myData_Gyro)

myData_Full <- left_join(myData_Acc, myData_Gyro, by = c("epoch", "user_id", "exp_id")) %>% 
        dplyr::select(-activity.y, -sample.y, -n.y) 

myData_Full$activity.x <-  as.integer(myData_Full$activity.x)


#SOME REMAINING ISSUES
#You should maybe think about the following issues

#What do you do with the unlabelled epochs? That is, the epochs that are marked with -???

#For the competition submision you should only provide predictions for the epochs defined 
#in the sample submission file???

#Should you remove the epochs that do not consist of 128 samples (i.e., n ≠ 128 in the above 
#data frame)????

##Additional features 
entropy  <- function(x, nbreaks = nclass.Sturges(x)) {
  r = range(x)
  x_binned = findInterval(x, seq(r[1], r[2], len= nbreaks))
  h = tabulate(x_binned, nbins = nbreaks) # fast histogram
  p = h/sum(h)
  -sum(p[p>0] * log(p[p>0]))
}



##MODELIING 

data_to_work_with <- myData_Full%>% filter(n.x >100) %>% select( -c(1:3, sample.x, n.x)) %>% drop_na()

model1 <- lm(activity.x~., data = dplyr::select(myData_Full, -epoch, -user_id, -exp_id,
                                                -sample.x, -n.x))
summary(model1)

car::vif(model1)









