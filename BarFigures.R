library(ggplot2)
# Read the data
Data = read.csv('AllFinalResultsOnset+5sec_200msec_bins.csv')

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
  geom_point(data=Data,aes(color=Trial),show.legend = FALSE)+
  geom_line(data=Data,aes(color=Trial,group=Trial),show.legend = FALSE) +
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
  geom_point(data=Data,aes(color=Trial),show.legend = FALSE)+
  geom_line(data=Data,aes(color=Trial,group=Trial),show.legend = FALSE) +
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
  geom_point(data=Data,aes(color=Trial),show.legend = FALSE)+
  geom_line(data=Data,aes(color=Trial,group=Trial),show.legend = FALSE) +
  labs(title="P-Value Metrics", x="Session Type", y = "Number of cells")+
  theme_classic()+
  scale_fill_manual(values=c('Red','Blue'))+
  theme(plot.title = element_text(hjust = 0.5))

