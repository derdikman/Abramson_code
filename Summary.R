
library(ggplot2)
# Read the data
Data = read.csv('AllFinalResultsOnset+5sec_200msec_bins.csv')

CellType.Trials = aggregate(ByCellType ~ Cell + SessionType + Trial, data = Data, sum) # sum all ByCellType per the Cell, trial and Session Type
CellType.Animal = aggregate(ByCellType ~ Cell + SessionType + Animal, data = Data, sum) # sum all ByCellType per the Cell, animal and Session Type
CellType.All.sum = aggregate(ByCellType ~ Cell + SessionType, data = Data, sum) # sum all ByCellType per the Cell, animal and Session Type
CellType.All.mean = aggregate(ByCellType ~ Cell + SessionType, data = Data, mean) # sum all ByCellType per the Cell, animal and Session Type

Fit.Trials = aggregate(ByFit ~ Cell + SessionType + Trial, data = Data, sum) # sum all ByFit per the Cell, trial and Session Type
Fit.Animal = aggregate(ByFit ~ Cell + SessionType + Animal, data = Data, sum) # sum all ByFit per the Cell, animal and Session Type
Fit.All.sum = aggregate(ByFit ~ Cell + SessionType, data = Data, sum) # sum all ByFit per the Cell, animal and Session Type
Fit.All.mean = aggregate(ByFit ~ Cell + SessionType, data = Data, mean) # sum all ByFit per the Cell, animal and Session Type

PValue.Trials = aggregate(ByPValue ~ Cell + SessionType + Trial, data = Data, sum) # sum all PValue per the Cell, trial and Session Type
PValue.Animal = aggregate(ByPValue ~ Cell + SessionType + Animal, data = Data, sum) # sum all PValue per the Cell, animal and Session Type
PValue.All.sum = aggregate(ByPValue ~ Cell + SessionType, data = Data, sum) # sum all PValue per the Cell, animal and Session Type
PValue.All.mean = aggregate(ByPValue ~ Cell + SessionType, data = Data, mean) # sum all PValue per the Cell, animal and Session Type

CellType.Totals = aggregate(ByCellType ~ SessionType + Animal, data = Data, sum) # sum # of ByCellType per the Cell and Session Type
Fit.Totals = aggregate(ByFit ~ SessionType + Animal, data = Data, sum) # sum # of ByFit per the Cell and Session Type
PValue.Totals = aggregate(ByPValue ~ SessionType + Animal, data = Data, sum) # sum # of ByPValue per the Cell and Session Type

with(Data,table(Animal,SessionType)) # Show number of experiments for every animal


#  CellType by sum of all session and animals
a.sum = with(CellType.All.sum,ByCellType[SessionType=="FixedDistance"]) # A vector with the cells type counting for FixedDistance
b.sum  = with(CellType.All.sum,ByCellType[SessionType=="FixedTime"]) # # A vector with the cells type counting for FixedTime
c.sum = as.table(rbind(a.sum,b.sum))
dimnames(c.sum) <- list(SessionType = c("Fixed Distance", "Fixed Time"), Cell = c("Distance", "Time"))

# CellType by mean of all session and animals
a.mean = with(CellType.All.mean,ByCellType[SessionType=="FixedDistance"]) # A vector with the cells type counting for FixedDistance
b.mean  = with(CellType.All.mean,ByCellType[SessionType=="FixedTime"]) # # A vector with the cells type counting for FixedTime
c.mean = as.table(rbind(a.mean,b.mean))
dimnames(c.mean) <- list(SessionType = c("Fixed Distance", "Fixed Time"), Cell = c("Distance", "Time"))

barplot(t(c.mean),legend=T,beside=T,main='CellType metrics (mean of all session and animals)')
# Pearson's Chi-squared Test for the mean count
chisq.test(c.mean)

barplot(t(c.sum),legend=T,beside=T,main='CellType metrics (sum of all session & animals)')
# Pearson's Chi-squared Test for the sum count
chisq.test(c.sum)


# CellType Per animal analysis
bk35=with(CellType.Animal,ByCellType[Animal=="bk35"])
bk35.c=as.table(rbind(bk35[1:2],bk35[3:4]))

bk41=with(CellType.Animal,ByCellType[Animal=="bk41"])
bk41.c=as.table(rbind(bk41[1:2],bk41[3:4]))

bk49=with(CellType.Animal,ByCellType[Animal=="bk49"])
bk49.c=as.table(rbind(bk49[1:2],bk49[3:4]))

bk45=with(CellType.Animal,ByCellType[Animal=="bk45"])
bk45.c=as.table(rbind(bk45[1:2]))

bk26=with(CellType.Animal,ByCellType[Animal=="bk26"])
bk26.c=as.table(rbind(bk26[1:2]))

bk33=with(CellType.Animal,ByCellType[Animal=="bk33"])
bk33.c=as.table(rbind(bk33[1:2]))

# Pearson's Chi-squared Test for Count Data, per animal
chisq.test(bk26.c)
chisq.test(bk33.c)
chisq.test(bk35.c)
chisq.test(bk41.c)
chisq.test(bk45.c)
chisq.test(bk49.c)

# Fit metrics by sum of all session and animals
a.sum = with(Fit.All.sum,ByFit[SessionType=="FixedDistance"]) # A vector with the cells type counting for FixedDistance
b.sum  = with(Fit.All.sum,ByFit[SessionType=="FixedTime"]) # # A vector with the cells type counting for FixedTime
c.sum = as.table(rbind(a.mean,b.mean))
dimnames(c.sum) <- list(SessionType = c("Fixed Distance", "Fixed Time"), Cell = c("Distance", "Time"))

# Fit metrics by mean of all session and animals
a.mean = with(Fit.All.mean,ByFit[SessionType=="FixedDistance"]) # A vector with the cells type counting for FixedDistance
b.mean  = with(Fit.All.mean,ByFit[SessionType=="FixedTime"]) # # A vector with the cells type counting for FixedTime
c.mean = as.table(rbind(a.mean,b.mean))
dimnames(c.mean) <- list(SessionType = c("Fixed Distance", "Fixed Time"), Cell = c("Distance", "Time"))

barplot(t(c.mean),legend=T,beside=T,main='Fit metrics (mean of all session and animals)')
# Pearson's Chi-squared Test for the mean count
chisq.test(c.mean)

barplot(t(c.sum),legend=T,beside=T,main='Fit metrics (sum of all session & animals)')
# Pearson's Chi-squared Test for the sum count
chisq.test(c.sum)

# PValue metrics by sum of all session and animals
a.sum = with(PValue.All.sum,ByFit[SessionType=="FixedDistance"]) # A vector with the cells type counting for FixedDistance
b.sum  = with(PValue.All.sum,ByFit[SessionType=="FixedTime"]) # # A vector with the cells type counting for FixedTime
c.sum = as.table(rbind(a.mean,b.mean))
dimnames(c.sum) <- list(SessionType = c("Fixed Distance", "Fixed Time"), Cell = c("Distance", "Time"))

# PValue metrics by mean of all session and animals  
a.mean = with(PValue.All.mean,ByPValue[SessionType=="FixedDistance"]) # A vector with the cells type counting for FixedDistance
b.mean  = with(PValue.All.mean,ByPValue[SessionType=="FixedTime"]) # # A vector with the cells type counting for FixedTime
c.mean = as.table(rbind(a.mean,b.mean))
dimnames(c.mean) <- list(SessionType = c("Fixed Distance", "Fixed Time"), Cell = c("Distance", "Time"))

barplot(t(c.mean),legend=T,beside=T,main='P Value metrics (mean of all session and animals)')
# Pearson's Chi-squared Test for the mean count
chisq.test(c.mean)

barplot(t(c.sum),legend=T,beside=T,main='P Value metrics (sum of all session & animals)')
# Pearson's Chi-squared Test for the sum count
chisq.test(c.sum)













