% Present the population analysis vs shuffle in voilin graphs
% The variable shown is the ratio between time and distance cells numbers on
% the time fixed trials and the ration between distance and time cells
% numbers on the fixed distance trials.
%
% Red line = median of the shuffled distribution 
% Black line = average of the shuffled distribution
% Blue line = actual ratio
%
% The AnalysisDuration='16sec' was chosen for the shuffles in order to
% provide the same duration for all trails (time trials duration=16sec,
% distance trial durations are longer).

clear;
close all;

AnalysisType = 'Onset'; % Onset based analysis
% AnalysisType = 'Peak'; % Peak based analysis
% AnalysisDuration = '+5sec_100msBins'
AnalysisDuration = '16sec';

% Load the full analysis r(1)=ByPValue, r(2)=ByFit r(3)=ByCellType
[r(1) r(2) r(3)]=AnalyzePopulation(['AnalyzedData_' AnalysisType '+5sec_100msBins'],0); % Population results 

% Load the truncated analysis to 16sec, t(1)=ByPValue, t(2)=ByFit t(3)=ByCellType
[t(1) t(2) t(3)]=AnalyzePopulation(['AnalyzedData_' AnalysisType '16sec'],0); % Population results at upto 16sec trials

% A vector of trials numbers the animal attended
bk49 = [8:9,17:18]; 
bk41 = [6:7,12:13];
bk35 = [5,10:11];
bk26 = [1:3];
bk33 = [4];
bk45 = [14:16];
all = [1:18];

n=1000; % Number of shuffles

for i=1:n
    [p{i}(1),p{i}(2),p{i}(3)] =AnalyzePopulation(['AnalyzedData_' AnalysisType AnalysisDuration],1,all); % Shuffle the trial type (shuffle=1 in the parameters)
end

% for i=1:n %Collect the cells encoding statistics
%     PValueT(i)   = p{i}(1).Ratio(1);
%     FitT(i)      = p{i}(2).Ratio(1);
%     CellTypeT(i) = p{i}(3).Ratio(1); 
%     
%     PValueD(i)   = 1/p{i}(1).Ratio(2);
%     FitD(i)      = 1/p{i}(2).Ratio(2);
%     CellTypeD(i) = 1/p{i}(3).Ratio(2); 
% end
for i=1:n %Collect the cells encoding statistics
    PValueT(i)   = p{i}(1).MajIndex(1);
    FitT(i)      = p{i}(2).MajIndex(1);
    CellTypeT(i) = p{i}(3).MajIndex(1); 
    
    PValueD(i)   = p{i}(1).MajIndex(2);
    FitD(i)      = p{i}(2).MajIndex(2);
    CellTypeD(i) = p{i}(3).MajIndex(2); 
end


%----------------------Histograms-----------------------------------------
h1=figure(1);
subplot(3,2,1);
h1=histogram(PValueT);
title('PValue Fixed Distance');
xlim([0,12]);
xticks([0,0.5,1:12]);
set(gca,'Xscale','log');
hold on;
plot(t(1).Ratio(1),0,'ok','markerfacecolor','b'); % 16sec Truncted results
plot(r(1).Ratio(1),0,'ok','markerfacecolor','r'); % Full trial results
[m ix] = max(h1.Values);
av1=(h1.BinEdges(ix)+h1.BinEdges(ix+1))/2;

subplot(3,2,2);
h2=histogram(PValueD);
title('PValue Fixed Time');
xlim([0,12]);
xticks([0,0.5,1:12]);
set(gca,'Xscale','log');
hold on;
plot(t(1).Ratio(2),0,'ok','markerfacecolor','b'); % 16sec Truncted results
plot(r(1).Ratio(2),0,'ok','markerfacecolor','r'); % Full trial results
[m ix] = max(h2.Values);
av2 =(h2.BinEdges(ix)+h2.BinEdges(ix+1))/2;

subplot(3,2,3);
h3=histogram(FitT);
title('Fit Fixed Distance');
xlim([0,4]);
xticks([0,0.5,1:4]);
set(gca,'Xscale','log');
hold on;
plot(t(2).Ratio(1),0,'ok','markerfacecolor','b'); % 16sec Truncted results
plot(r(2).Ratio(1),0,'ok','markerfacecolor','r'); % Full trail results
[m ix] = max(h3.Values);
av3=(h3.BinEdges(ix)+h3.BinEdges(ix+1))/2;

subplot(3,2,4);
h4=histogram(FitD);
title('Fit Fixed Time');
xlim([0,4]);
xticks([0,0.5,1:4]);
set(gca,'Xscale','log');
hold on;
plot(t(2).Ratio(2),0,'ok','markerfacecolor','b'); % 16sec Truncted results
plot(r(2).Ratio(2),0,'ok','markerfacecolor','r'); % Full trail results
[m ix] = max(h4.Values);
av4=(h4.BinEdges(ix)+h4.BinEdges(ix+1))/2;

subplot(3,2,5);
h5=histogram(CellTypeT);
title('CellType Fixed Distance');
xlim([0,3]);
xticks([0,0.5,1:3]);
set(gca,'Xscale','log');
hold on;
plot(t(3).Ratio(1),0,'ok','markerfacecolor','b'); % 16sec Truncted results
plot(r(3).Ratio(1),0,'ok','markerfacecolor','r'); % Full trail results
[m ix] = max(h5.Values);
av5=(h5.BinEdges(ix)+h5.BinEdges(ix+1))/2;

subplot(3,2,6);
h6=histogram(CellTypeD);
title('CellType Fixed Time');
xlim([0.05,3]);
xticks([0.05,0.5,1:3]);
set(gca,'Xscale','log');
hold on;
plot(t(3).Ratio(2),0,'ok','markerfacecolor','b'); % 16sec Truncted results
plot(r(3).Ratio(2),0,'ok','markerfacecolor','r'); % Full trail results
[m ix] = max(h6.Values);
av6=(h6.BinEdges(ix)+h6.BinEdges(ix+1))/2;

sgtitle(['Distribution of ',num2str(n),' Trial Type shuffles of 16sec truncated runs, Peak analysis']);


%----------- Violin graphes -----------------------------------------
h2=figure(2);
h2.Position = [2340 411 1242 420];

subplot(1,3,1);
violin([PValueT',PValueD'],'xlabel',{'Fixed-Distance','Fixed-Time'});
ylim([-1,1]);
yticks([-1,0,1]);
% Present the actual ratio of the 16sec truncated analysis
plot([0.75,1.25],[t(1).MajIndex(1),t(1).MajIndex(1)],'-b');
plot([1.75,2.25],[t(1).MajIndex(2),t(1).MajIndex(2)],'-b');

% Present the actual ratio of the full analysis
% plot([0.75,1.25],[r(1).MajIndex(1),r(1).MajIndex(1)],'-g');
% plot([1.75,2.25],[r(1).MajIndex(2),r(1).MajIndex(2)],'-g');
legend('off');
title('PValue');

subplot(1,3,2);
violin([FitT',FitD'],'xlabel',{'Fixed-Distance','Fixed-Time'});
ylim([-0.5,0.5]);
yticks([-0.5,0,0.5]);

% Present the actual ratio of the 16sec truncated analysis
plot([0.75,1.25],[t(2).MajIndex(1),t(2).MajIndex(1)],'-b');
plot([1.75,2.25],[t(2).MajIndex(2),t(2).MajIndex(2)],'-b');

% Present the actual ratio of the full analysis
% plot([0.75,1.25],[r(2).Ratio(1),r(2).Ratio(1)],'-g');
% plot([1.75,2.25],[1/r(2).Ratio(2),1/r(2).Ratio(2)],'-g');

legend('off');
title('Fit')


subplot(1,3,3);
violin([CellTypeT',CellTypeD'],'xlabel',{'Fixed-Distance','Fixed-Time'});
ylim([-0.3,0.3]);
yticks([-0.3,0,0.3]);
% Present the actual ratio of the 16sec truncated analysis
plot([0.75,1.25],[t(3).MajIndex(1),t(3).MajIndex(1)],'-b');
plot([1.75,2.25],[t(3).MajIndex(2),t(3).MajIndex(2)],'-b');

% Present the actual ratio of the full analysis
% plot([0.75,1.25],[r(3).Ratio(1),r(3).Ratio(1)],'-g');
% plot([1.75,2.25],[1/r(3).Ratio(2),1/r(3).Ratio(2)],'-g');

legend('off');
title('CellType');

sgtitle([num2str(n),' Trial Type shuffles of 16sec truncated runs, ' AnalysisType ' analysis']);


disp([av1 av2 av1*av2]);
disp([av3 av4 av3*av4]);
disp([av5 av6 av5*av6]);

disp('Full track results, by CellType :')
disp(['              FDist FTime']);
disp(sprintf('%s  %d  %d','DistanceCells ',r(3).Distance));
disp(sprintf('%s  %d  %d','TimeCells     ',r(3).Time));
disp(sprintf('%s  %d','Total Cells :   ' ,sum(r(3).Time)+sum(r(3).Distance)));
disp('   ');
disp('16sec truncated track results, by CellType :')
disp(['              FDist FTime']);
disp(sprintf('%s  %d  %d','DistanceCells ',t(3).Distance));
disp(sprintf('%s  %d  %d','TimeCells     ',t(3).Time));
disp(sprintf('%s  %d','Total Cells :   ' ,sum(t(3).Time)+sum(t(3).Distance)));



