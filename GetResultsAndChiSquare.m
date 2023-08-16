% Calculate the Chi-Square results for each animal and for all animals


Animals = {'BK49','BK41','BK35','BK45','BK33','BK26','All'};
Trials = { [8,9,17,18], [6,7,12,13],[5,10,11],[14,15,16], [4],[1,2,3],[1:18]}; % The enumerator of the respective trial

file = uigetfile;
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation(file,0);

for Animal=1:7
    
    CellType.Time = [sum(ByCellType.TimeCells(TrialType<10 & ismember(TrialType,Trials{Animal}))),sum(ByCellType.TimeCells(TrialType>9 & ismember(TrialType,Trials{Animal})))];
    CellType.Distance = [sum(ByCellType.DistanceCells(TrialType<10 & ismember(TrialType,Trials{Animal}))),sum(ByCellType.DistanceCells(TrialType>9 & ismember(TrialType,Trials{Animal})))];
    
    Fit.Time = [sum(ByFit.TimeCells(TrialType<10 & ismember(TrialType,Trials{Animal}))),sum(ByFit.TimeCells(TrialType>9 & ismember(TrialType,Trials{Animal})))];
    Fit.Distance = [sum(ByFit.DistanceCells(TrialType<10 & ismember(TrialType,Trials{Animal}))),sum(ByFit.DistanceCells(TrialType>9 & ismember(TrialType,Trials{Animal})))];
    
    PValue.Time = [sum(ByPValue.TimeCells(TrialType<10 & ismember(TrialType,Trials{Animal}))),sum(ByPValue.TimeCells(TrialType>9 & ismember(TrialType,Trials{Animal})))];
    PValue.Distance = [sum(ByPValue.DistanceCells(TrialType<10 & ismember(TrialType,Trials{Animal}))),sum(ByPValue.DistanceCells(TrialType>9 & ismember(TrialType,Trials{Animal})))];
    
    [pval(1), x2a(1)] = chisquarecont([CellType.Time;CellType.Distance]);
    [pval(2), x2a(2)] = chisquarecont([Fit.Time(1:2);Fit.Distance(1:2)]);
    [pval(3), x2a(3)] = chisquarecont([PValue.Time;PValue.Distance]);
    
     
    tCellType=table(CellType.Distance',CellType.Time','VariableNames', {'FD','FT'} );
    disp(Animals(Animal));
    if ~isnan(pval(1))
        disp('       P    Chi-Square');
        disp([pval(1), x2a(1)]);
        disp([pval(2), x2a(2)]);
        disp([pval(3), x2a(3)]);
        disp(tCellType);
    end
    
    
    
end
