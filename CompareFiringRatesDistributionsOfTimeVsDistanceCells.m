% Compares the firing rates distributions of the Time vs Distance cells


% Load the desired mat file containing R1-R18 results summary
file = 'AnalyzedData_Onset+5sec_100msBins.mat';
load(file);
[ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation(file,0,[1:18]); 

Peaks=max(R(:,30:32)')'; % Peak firing rates for all cells
Avg=R(:,end);

figure(1);
clf;
subplot(1,2,1);

T = Peaks(ByCellType.TimeCells);
D = Peaks(ByCellType.DistanceCells);

[h1,edges] = histcounts(T, 10);
[h2,edges] = histcounts(D, 10);
ctrs = edges(1)+(1:length(edges)-1).*diff(edges);   % Create Centres
bar(ctrs, [h1/sum(h1)*100 ;h2/sum(h2)*100]')


xlabel('Hz');
ylabel('% of cells');
legend('Time cells', 'Distance cells');

% Perform the Kolmogorov-Smirnov test
[h, pValue, ksStat] = kstest2(h1/sum(h1), h2/sum(h2));

title('Max firing rates of Time and Distance cells');
subtitle(['p-value (KS test): ' num2str(pValue,3)]);

subplot(1,2,2);

T = Avg(ByCellType.TimeCells);
D = Avg(ByCellType.DistanceCells);

[h1,edges] = histcounts(T, 10);
[h2,edges] = histcounts(D, 10);
ctrs = edges(1)+(1:length(edges)-1).*diff(edges);   % Create Centres
bar(ctrs, [h1/sum(h1)*100 ;h2/sum(h2)*100]')


xlabel('Hz');
ylabel('% of cells');
legend('Time cells', 'Distance cells');

% Perform the Kolmogorov-Smirnov test
[h, pValue, ksStat] = kstest2(h1/sum(h1), h2/sum(h2));

title('Average firing rates of Time and Distance cells');
subtitle(['p-value (KS test): ' num2str(pValue,3)]);

