KrausData :  .mat files containing the data from each of the trials

AnalyzeAllData.m : matlab function running analysis of all session and creates 'AnalyzedData.mat'
RunAllSessionsTimeVsDistance.m : matlab script running the analysis of all trials (called from AnalyzedAllData.m)
NeuronFiringTimeVsDistance.m : matlab script analyzing a given trial (called from TunAllSessionsTimeVsDistance.m)

GetResultAndChiSquare.m :  Calculate the Chi-Square results for each animal and for all animals
AnalyzePopulation.m : Based on the analysis file given, calculate the number of cells in every category
CompareFiringRatesDistributionsOFTimeVsDistance.m : Compares the firing rates distributions of the Time vs Distance cells
CompareMaxFiringRatesVsVelocityPerCellType.m : The script plots the averaged max firing rate of each run vs the velocity for time and distance cells.
CompareMaxFiringRatesVsVelocityPerTrialType.m : The script plots the averaged max firing rate of each run vs the velocity in Fixed-Time and Fixed-Distance experiments.
ROC.m : Plot the Receiver Operating Characteristic and Youden index 
ShuffleHistogram :  Present the population analysis vs shuffle in voilin graphs
 
Results mat files : 
AllFinalResultsOnset+5sec_100mSec_bins.mat : Analysis results with 100mSec bins, Onset analysis
AnalyzedData_Onset16sec.mat : Analysis results with 100mSec bins, Onset analysis, truncated to 16sec
AllFinalResultsPeak+5sec_100mSec_bins.mat : Analysis results with 100mSec bins, Peak analysis
AnalyzedData_Peak16sec.mat : Analysis results with 100mSec bins, Peak analysis, truncated to 16sec

BarFigures.R : an R script using the above csv file and generating the paper bar figures
