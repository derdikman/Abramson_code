library(ggplot2)
library(tidyverse)  
# Read the data
Data = read.csv('G:/Dropbox (Personal)/Brain/Project/AllDataOnset+5with200mSecBins.csv')

df <- Data
df$SessionType <- as.factor(df$SessionType)
df$Cell <- as.factor(df$Valid)
head(df)

df2 = filter(df,Valid == TRUE & (Animal=='bk35' |Animal=='bk41' | Animal=='bk49'))
#df2 = filter(df,df$Valid == TRUE, Animal == 'bk35')

# Plot the Mean Time for ByCellType classification at FTime and FDistance experiments
df2 %>%
  ggplot(aes(x=Mtime,fill=ByCellType)) +
  facet_wrap(~ByCellType+~SessionType, strip.position = "bottom") +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity',bins=10) +
  scale_fill_manual(values=c("Red", "Blue"))
#  labs(fill="")

# Plot the Mean Distance for ByCellType classification at FTime and FDistance experiments
df2 %>%
  ggplot(aes(x=MDist,fill=ByCellType)) +
  facet_wrap(~ByCellType+~SessionType, strip.position = "bottom") +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity',bins=10) +
  scale_fill_manual(values=c("Red", "Blue"))
#  labs(fill="")

###########################################################
###########################################################
###########################################################
# Plot the Mean Time for all cells and all animals at FTime and FDistance experiments
#df3 = filter(df,df$Valid == TRUE & (df$Animal=='bk35' |df$Animal=='bk41' | df$Animal=='bk49'))
#df3 = filter(df,df$Valid == TRUE & (df$Animal=='bk26' |df$Animal=='bk33' | df$Animal=='bk45'))
df3 = filter(df,df$Valid == TRUE)


dfFT = filter(df3, df3$Valid == TRUE & df3$SessionType=="FixedTime")
dfFD = filter(df3,df3$Valid == TRUE & df3$SessionType=="FixedDistance")
mean(dfFD$MDist)-mean(dfFT$MDist)
median(dfFD$MDist)-median(dfFT$MDist)

df3 %>%
  ggplot(aes(x=Mtime,fill=SessionType)) +
  facet_wrap(~SessionType+Animal, strip.position = "bottom") +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity',bins=10) +
  scale_fill_manual(values=c("Red", "Blue"))
#  labs(fill="")

df3 %>%
  ggplot( aes(x=SessionType, y=Mtime, fill=SessionType)) +
  facet_wrap(~SessionType+Animal, strip.position = "bottom") +
    geom_boxplot() +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  # theme_ipsum() +
  theme(
    plot.title = element_text(size=11)
  ) +
  scale_fill_manual(values=c("Red", "Blue")) +
  ggtitle("A boxplot with jitter") +
  xlab("")

# Plot the Mean Distance for all cells and all animals at FTime and FDistance experiments
df3 %>%
  ggplot(aes(x=MDist,fill=SessionType)) +
  facet_wrap(~Animal, strip.position = "bottom",ncol=2) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity',bins=20) +
  scale_fill_manual(values=c("Red", "Blue"))


df3 %>%
  ggplot( aes(x=SessionType, y=MDist, fill=SessionType)) +
  facet_wrap(~Animal, strip.position = "bottom",ncol = 3) +
  geom_boxplot() +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  # theme_ipsum() +
  theme(
    plot.title = element_text(size=11)
  ) +
  scale_fill_manual(values=c("Red", "Blue")) +
  ggtitle("A boxplot with jitter") +
  xlab("")


###########################################################
###########################################################

# Plot
df2 %>%
  ggplot( aes(x=SessionType, y=Mtime, fill=ByCellType)) +
  geom_boxplot() +
 # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
 # theme_ipsum() +
  theme(
    plot.title = element_text(size=11)
  ) +
  scale_fill_manual(values=c("Red", "Blue")) +
  ggtitle("A boxplot with jitter") +
  xlab("")

# Aggregate all experiments CellType-metrics at each session type and cell type , calculating the mean,min,max and sem
myData <- aggregate(Data$ByCellType,
                    by = list(Cell = Data$Cell, SessionType = Data$SessionType),
                    FUN = function(x) c(mean = mean(x), max=max(x), min=min(x), sem = sd(x)/sqrt(length(x)), n = length(x))
)

# Convert to a matrix
myData <- do.call(data.frame, myData)

Data$x.mean=Data$ByCellType # Added so we can use Data for the geom_point and geom_line in order to display all trials points

ggplot(myData,aes(x=Cell, y=x.mean, fill=Cell))+
  facet_wrap(~SessionType, strip.position = "bottom") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=x.mean+x.sem, ymax=x.mean-x.sem), width=.2, position=position_dodge(.9),size =0.9) +
  geom_point(data=Data,aes(color=Animal),show.legend = FALSE)+
  geom_line(data=Data,aes(color=Animal,group=Trial),show.legend = TRUE) +
  labs(title="CellType Metrics", x="Session Type", y = "Number of cells")+
  theme_classic()+
  scale_fill_manual(values=c('Red','Blue'))+
  theme(plot.title = element_text(hjust = 0.5))



# Aggregate all experiments Fit-metrics at each session type and cell type , calculating the mean,min,max and sem
myData <- aggregate(Data$ByFit,
                    by = list(Cell = Data$Cell, SessionType = Data$SessionType),
                    FUN = function(x) c(mean = mean(x), max=max(x), min=min(x), sem = sd(x)/sqrt(length(x)), n = length(x))
)

# Convert to a matrix
myData <- do.call(data.frame, myData)

Data$x.mean=Data$ByFit

ggplot(myData,aes(x=Cell, y=x.mean, fill=Cell))+
  facet_wrap(~SessionType, strip.position = "bottom") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=x.mean+x.sem, ymax=x.mean-x.sem), width=.2, position=position_dodge(.9),size =0.9) +
  geom_point(data=Data,aes(color=Animal),show.legend = FALSE)+
  geom_line(data=Data,aes(color=Animal,group=Trial),show.legend = TRUE) +
  labs(title="Fit Metrics", x="Session Type", y = "Number of cells")+
  theme_classic()+
  scale_fill_manual(values=c('Red','Blue'))+
  theme(plot.title = element_text(hjust = 0.5))


# Aggregate all experiments P-Value-metrics at each session type and cell type , calculating the mean,min,max and sem
myData <- aggregate(Data$ByPValue,
                    by = list(Cell = Data$Cell, SessionType = Data$SessionType),
                    FUN = function(x) c(mean = mean(x), max=max(x), min=min(x), sem = sd(x)/sqrt(length(x)), n = length(x))
)

# Convert to a matrix
myData <- do.call(data.frame, myData)

Data$x.mean=Data$ByPValue

ggplot(myData,aes(x=Cell, y=x.mean, fill=Cell))+
  facet_wrap(~SessionType, strip.position = "bottom") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=x.mean+x.sem, ymax=x.mean-x.sem), width=.2, position=position_dodge(.9),size =0.9) +
  geom_point(data=Data,aes(color=Animal),show.legend = FALSE)+
  geom_line(data=Data,aes(color=Animal,group=Trial),show.legend = TRUE) +
  labs(title="P-Value Metrics", x="Session Type", y = "Number of cells")+
  theme_classic()+
  scale_fill_manual(values=c('Red','Blue'))+
  theme(plot.title = element_text(hjust = 0.5))





