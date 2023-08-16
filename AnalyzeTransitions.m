
% Analyze the animals with transitions from fixed-time to fixed-distance
% and vice versa, to trace changes in encoding

% BK41
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsFirst20Runs.mat',0);
BK41(1,1)=(sum(ByCellType.DistanceCells(TrialType==6))/sum(ByCellType.TimeCells(TrialType==6)));
disp([sum(ByCellType.DistanceCells(TrialType==6)),sum(ByCellType.TimeCells(TrialType==6))]);
BK41(2,1)=(sum(ByCellType.TimeCells(TrialType==12))/sum(ByCellType.DistanceCells(TrialType==12)));
disp([sum(ByCellType.DistanceCells(TrialType==12)),sum(ByCellType.TimeCells(TrialType==12))]);
BK41(3,1)=(sum(ByCellType.TimeCells(TrialType==13))/sum(ByCellType.DistanceCells(TrialType==13)));
disp([sum(ByCellType.DistanceCells(TrialType==13)),sum(ByCellType.TimeCells(TrialType==13))]);
BK41(4,1)=(sum(ByCellType.DistanceCells(TrialType==7))/sum(ByCellType.TimeCells(TrialType==7)));
disp([sum(ByCellType.DistanceCells(TrialType==7)),sum(ByCellType.TimeCells(TrialType==7))]);
disp('------------------');
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsLast20Runs.mat',0);
BK41(1,2)=(sum(ByCellType.DistanceCells(TrialType==6))/sum(ByCellType.TimeCells(TrialType==6)));
disp([sum(ByCellType.DistanceCells(TrialType==6)),sum(ByCellType.TimeCells(TrialType==6))]);
BK41(2,2)=(sum(ByCellType.TimeCells(TrialType==12))/sum(ByCellType.DistanceCells(TrialType==12)));
disp([sum(ByCellType.DistanceCells(TrialType==12)),sum(ByCellType.TimeCells(TrialType==12))]);
BK41(3,2)=(sum(ByCellType.TimeCells(TrialType==13))/sum(ByCellType.DistanceCells(TrialType==13)));
disp([sum(ByCellType.DistanceCells(TrialType==13)),sum(ByCellType.TimeCells(TrialType==13))]);
BK41(4,2)=(sum(ByCellType.DistanceCells(TrialType==7))/sum(ByCellType.TimeCells(TrialType==7)));
disp([sum(ByCellType.DistanceCells(TrialType==7)),sum(ByCellType.TimeCells(TrialType==7))]);
disp('D-T-T-D');
disp(BK41);

% BK49
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsFirst20Runs.mat',0);
BK49(1,1)=(sum(ByCellType.DistanceCells(TrialType==8))/sum(ByCellType.TimeCells(TrialType==8)));
disp([sum(ByCellType.DistanceCells(TrialType==8)),sum(ByCellType.TimeCells(TrialType==8))]);
BK49(2,1)=(sum(ByCellType.DistanceCells(TrialType==9))/sum(ByCellType.TimeCells(TrialType==9)));
disp([sum(ByCellType.DistanceCells(TrialType==9)),sum(ByCellType.TimeCells(TrialType==9))]);
BK49(3,1)=(sum(ByCellType.TimeCells(TrialType==17))/sum(ByCellType.DistanceCells(TrialType==17)));
disp([sum(ByCellType.DistanceCells(TrialType==17)),sum(ByCellType.TimeCells(TrialType==17))]);
BK49(4,1)=(sum(ByCellType.TimeCells(TrialType==18))/sum(ByCellType.DistanceCells(TrialType==18)));
disp([sum(ByCellType.DistanceCells(TrialType==18)),sum(ByCellType.TimeCells(TrialType==18))]);
disp('------------------');
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsLast20Runs.mat',0);
BK49(1,2)=(sum(ByCellType.DistanceCells(TrialType==8))/sum(ByCellType.TimeCells(TrialType==8)));
disp([sum(ByCellType.DistanceCells(TrialType==8)),sum(ByCellType.TimeCells(TrialType==8))]);
BK49(2,2)=(sum(ByCellType.DistanceCells(TrialType==9))/sum(ByCellType.TimeCells(TrialType==9)));
disp([sum(ByCellType.DistanceCells(TrialType==9)),sum(ByCellType.TimeCells(TrialType==9))]);
BK49(3,2)=(sum(ByCellType.TimeCells(TrialType==17))/sum(ByCellType.DistanceCells(TrialType==17)));
disp([sum(ByCellType.DistanceCells(TrialType==17)),sum(ByCellType.TimeCells(TrialType==17))]);
BK49(4,2)=(sum(ByCellType.TimeCells(TrialType==18))/sum(ByCellType.DistanceCells(TrialType==18)));
disp([sum(ByCellType.DistanceCells(TrialType==18)),sum(ByCellType.TimeCells(TrialType==18))]);


disp('D-D-T-T');
disp(BK49);

% BK35
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsFirst20Runs.mat',0);
BK35(1,1)=(sum(ByCellType.DistanceCells(TrialType==5))/sum(ByCellType.TimeCells(TrialType==5)));
disp([sum(ByCellType.DistanceCells(TrialType==5)),sum(ByCellType.TimeCells(TrialType==5))]);
BK35(2,1)=(sum(ByCellType.TimeCells(TrialType==10))/sum(ByCellType.DistanceCells(TrialType==10)));
disp([sum(ByCellType.DistanceCells(TrialType==10)),sum(ByCellType.TimeCells(TrialType==10))]);
BK35(3,1)=(sum(ByCellType.TimeCells(TrialType==11))/sum(ByCellType.DistanceCells(TrialType==11)));
disp([sum(ByCellType.DistanceCells(TrialType==11)),sum(ByCellType.TimeCells(TrialType==11))]);
disp('------------------');
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsLast20Runs.mat',0);
BK35(1,2)=(sum(ByCellType.DistanceCells(TrialType==5))/sum(ByCellType.TimeCells(TrialType==5)));
disp([sum(ByCellType.DistanceCells(TrialType==5)),sum(ByCellType.TimeCells(TrialType==5))]);
BK35(2,2)=(sum(ByCellType.TimeCells(TrialType==10))/sum(ByCellType.DistanceCells(TrialType==10)));
disp([sum(ByCellType.DistanceCells(TrialType==10)),sum(ByCellType.TimeCells(TrialType==10))]);
BK35(3,2)=(sum(ByCellType.TimeCells(TrialType==11))/sum(ByCellType.DistanceCells(TrialType==11)));
disp([sum(ByCellType.DistanceCells(TrialType==11)),sum(ByCellType.TimeCells(TrialType==11))]);
disp('D-T-T');
disp(BK35)















% 
% [ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsLast20Runs.mat',0);
% disp(sum(ByCellType.DistanceCells(TrialType==6))/sum(ByCellType.TimeCells(TrialType==6)));
% disp(sum(ByCellType.TimeCells(TrialType==12))/sum(ByCellType.DistanceCells(TrialType==12)));
% disp(sum(ByCellType.TimeCells(TrialType==13))/sum(ByCellType.DistanceCells(TrialType==13)));
% disp(sum(ByCellType.DistanceCells(TrialType==7))/sum(ByCellType.TimeCells(TrialType==7)));
% 
% disp(sum(ByCellType.DistanceCells(TrialType==7))/sum(ByCellType.TimeCells(TrialType==7)));
% disp(sum(ByCellType.DistanceCells(TrialType==8))/sum(ByCellType.TimeCells(TrialType==8)));
% % Fixed distance to Fixed time
% disp(sum(ByCellType.TimeCells(TrialType==10))/sum(ByCellType.DistanceCells(TrialType==10)));
% disp(sum(ByCellType.TimeCells(TrialType==12))/sum(ByCellType.DistanceCells(TrialType==12)));
% 
% [ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBinsLast20Runs.mat',0);
% % Fixed time to Fixed distance
% disp(sum(ByCellType.DistanceCells(TrialType==7))/sum(ByCellType.TimeCells(TrialType==7)));
% disp(sum(ByCellType.DistanceCells(TrialType==8))/sum(ByCellType.TimeCells(TrialType==8)));
% % Fixed distance to Fixed time
% disp(sum(ByCellType.TimeCells(TrialType==10))/sum(ByCellType.DistanceCells(TrialType==10)));
% disp(sum(ByCellType.TimeCells(TrialType==12))/sum(ByCellType.DistanceCells(TrialType==12)));
% 
% [ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation('AnalyzedData_Onset+5sec_250msBins.mat',0);
% % Fixed time to Fixed distance
% disp(sum(ByCellType.DistanceCells(TrialType==7))/sum(ByCellType.TimeCells(TrialType==7)));
% disp(sum(ByCellType.DistanceCells(TrialType==8))/sum(ByCellType.TimeCells(TrialType==8)));
% % Fixed distance to Fixed time
% disp(sum(ByCellType.TimeCells(TrialType==10))/sum(ByCellType.DistanceCells(TrialType==10)));
% disp(sum(ByCellType.TimeCells(TrialType==12))/sum(ByCellType.DistanceCells(TrialType==12)));