% Plot the Receiver Operating Characteristic and Youden index 


file = uigetfile;
load(file);

minFiringRate=0.5; % *hz) Minimum firing rate to consider
maxFiringRate=15; % (hz) Maximum firing rate to consider
MinTime = 1; % Do not analyze cells with peak at time less than MinTime from treadmill start (Sec)
MinValidRunsPercentage = 0.1; % minimum required runs for analysis

TrialType = R(:,1); % 1-9 are fixed-Distance trials, 10-18 are fixed-Time trials
ValidRuns=R(:,4);
MeanTime=R(:,10);
MaxFR = max(R(1:end,28:30)'); % Maxmimum firing rate of the 3 velocity groups

% Validity criteria 
Valid = (ValidRuns>MinValidRunsPercentage) & (MaxFR>minFiringRate)' & (MaxFR<maxFiringRate)' & MeanTime>MinTime;
FixedDistance=TrialType<10; % The fixed distance trials
FixedTime=TrialType>=10; % The fixed time trials

% CellType metrics is based on the variance from the linear fit 
CellType = (R(:,20)-R(:,24))./(R(:,20)+R(:,24));
cnt = 0;
TPR=[];
FPR=[];
TNR=[];
FNR=[];
for CellTypeThreshold=1:-0.01:-1
    TimeCells= (CellType<CellTypeThreshold & Valid);
    DistanceCells= (CellType<CellTypeThreshold & Valid);
    Time=    [sum(CellType>=CellTypeThreshold  & FixedDistance & Valid),sum(CellType>=CellTypeThreshold  & FixedTime & Valid)];
    Distance=[sum(CellType<CellTypeThreshold & FixedDistance & Valid),sum(CellType<CellTypeThreshold & FixedTime & Valid)];
    cnt = cnt +1;
    Th(cnt)=CellTypeThreshold;
    TPR(cnt) = Distance(1)/(Distance(1)+Time(1)); % cells classified as distance in the fixed-distance trials
    TNR(cnt) = 1-TPR(cnt);
    FPR(cnt) = Distance(2)/(Distance(2)+Time(2)); % cells classified as distance in the fixed-time trials
    FNR(cnt) = 1-FPR(cnt);
    Youden(cnt) = TPR(cnt)-FPR(cnt);%Distance(1)/(Distance(1)+Time(2))+Time(1)/(Time(1)+Distance(2))-1;
%     disp([CellTypeThreshold TPR(cnt) FPR(cnt) Youden(cnt)]);
end
fig=figure(1);
clf;
fig.Position=[2056 241 601 463];
p1=plot(FPR,TPR,'-r');
hold on;
p2=plot(FPR(101),TPR(101),'ok','MarkerFaceColor', 'k','DisplayName', '-0.01');
p3=plot(FPR(102),TPR(102),'ok','MarkerFaceColor', 'g','DisplayName', '0');
p4=plot(FPR(100),TPR(100),'ok','MarkerFaceColor', 'b','DisplayName', '0.01');
xticks([0 0.5 1]);
yticks([0 0.5 1]);
xlabel('False Positive Rate','FontSize',14,'FontWeight', 'bold');
ylabel('True Positive Rate','FontSize',14,'FontWeight', 'bold');
legend([p4,p3,p2],'Location', 'southeast');
title('ROC curve for Distance Session predictor','FontWeight', 'bold','FontSize',16 ); 

figure(2);
plot([1:-0.01:-1],Youden,'-*');
title('Youden index'); 
xlabel('Threshold');
ylabel('Index');




