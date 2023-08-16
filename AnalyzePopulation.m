% Based on the analysis file given, calculate the number of cells in every category

function  [ByPValue,ByFit,ByCellType,TrialType]= AnalyzePopulation(file,Shuffle, Trials)
load(file);
if nargin<2
    Shuffle=0;
end
if nargin<3
    Trials=[1:18];
end

ShowBar = 0; % 0= do not show, 1=show the bar graph
minFiringRate=0.5; % *hz) Minimum firing rate to consider
maxFiringRate=15; % (hz) Maximum firing rate to consider
MinTime = 1; % Do not analyze cells with peak at time less than MinTime from treadmill start (Sec)
MinValidRunsPercentage = 0.1; % minimum required runs for analysis
FitThreshold=[0.5 100]; % Thresholds for the Fit metrics
CellTypeThreshold=0.0; % Thresholds for the CellType metrics
PvalueThreshold=0.05; % Thresholds for the PValue metrics

TrialType = R(:,1); % 1-9 are fixed-Distance trials, 10-18 are fixed-Time trials
ValidRuns=R(:,4);
MeanTime=R(:,10);
MaxFR = max(R(1:end,28:30)'); % Maxmimum firing rate of the 3 velocity groups

% Validity criteria 
Valid = (ValidRuns>MinValidRunsPercentage) & (MaxFR>minFiringRate)' & (MaxFR<maxFiringRate)' & MeanTime>MinTime & ismember(TrialType,Trials);
FixedDistance=TrialType<10; % The fixed distance trials
FixedTime=TrialType>=10; % The fixed time trials

% If shuffle, then randomize the Trial type (either Distance or Time)
if Shuffle==1
    NFD = sum(FixedDistance & Valid); % Number of valid cells in the Fixed Distance trials
    NFT = sum(FixedTime & Valid); % Number of valid cells in the Fixed Time trials
    FixedDistance(1:end)=0;
    FixedTime(1:end)=0;
    TrialType(1:end)=~Valid; %Use the TrialType vector to indicate which cells were chosen (0=free)
    for j=1:NFD % Select the cells for the (pretended) Fixed-Distance trials 
        ix = randi(length(Valid));
        while TrialType(ix)>0
            ix = randi(length(Valid));
        end
        TrialType(ix) = 9; % Fixed Distance indicator
    end
    for j=1:NFT % Select the cells for the (pretended) Fixed-Distance trials 
        ix = randi(length(Valid));
        while TrialType(ix)>0
            ix = randi(length(Valid));
        end
        TrialType(ix) = 18; % Fixed Time indicator
    end
    FixedDistance=TrialType==9;
    FixedTime=TrialType==18;
end

% P-value metrics is the p-value of rejecting the null hypothesis that
% there is no linear relation between distance and velcotiry or time and
% inverse velocity.
Ps_v = R(:,16);
Pt_rv = R(:,22);
ByPValue.TimeCells = (Ps_v<PvalueThreshold & Pt_rv>=PvalueThreshold & Valid);
ByPValue.DistanceCells = (Ps_v>=PvalueThreshold & Pt_rv<PvalueThreshold & Valid);
ByPValue.Time =[sum(Ps_v<PvalueThreshold & Pt_rv>=PvalueThreshold & FixedDistance & Valid),sum(Ps_v<PvalueThreshold  & Pt_rv>=PvalueThreshold & FixedTime & Valid)];
ByPValue.Distance=[sum(Ps_v>=PvalueThreshold & Pt_rv<PvalueThreshold & FixedDistance & Valid),sum(Ps_v>=PvalueThreshold & Pt_rv<PvalueThreshold & FixedTime & Valid)];
ByPValue.Total = [sum(ByPValue.Time),sum(ByPValue.Distance)];
% Ratio between distance and time cells in the distance trials and between time and distance in the time trials.
% If the ratio >1 then there is a relation between the cells' encoding and the trial type (either fixed distance or fixed time)
% The Majority Index is -1 to 1 (-1 if all cells are Time cells, 1 if all cells are distance cells)
ByPValue.Ratio = [ByPValue.Distance(1)/ByPValue.Time(1),ByPValue.Time(2)/ByPValue.Distance(2)]; 
ByPValue.MajIndex = [(ByPValue.Distance(1)-ByPValue.Time(1))/(ByPValue.Distance(1)+ByPValue.Time(1)),(ByPValue.Distance(2)-ByPValue.Time(2))/(ByPValue.Distance(2)+ByPValue.Time(2))]; 

% ByPValue.Both=[sum(Ps_v<0.05 & Pt_rv<0.05 & FixedDistance & Valid),sum(Ps_v<0.05 & Pt_rv<0.05 & FixedTime & Valid)];

% Fit metrics is based on the ratio between the linear fit slope of distance vs Velocity 
% or Time vs inverse velocity, compared to the mean distance or time
FitT=2*R(:,7)./R(:,10)-1;
FitD=1-2*R(:,8)./R(:,12);
ByFit.TimeCells= (FitT>FitThreshold(1) & FitT<FitThreshold(2) & Valid);
ByFit.DistanceCells= (FitD<-FitThreshold(1) & FitD>-FitThreshold(2) & Valid);
ByFit.Time=    [sum(FitT>FitThreshold(1) & FitT<FitThreshold(2)   & FixedDistance & Valid),sum(FitT>FitThreshold  & FitT<FitThreshold(2) & FixedTime & Valid)];
ByFit.Distance=[sum(FitD<-FitThreshold(1) & FitD>-FitThreshold(2) & FixedDistance & Valid),sum(FitD<-FitThreshold & FitD>-FitThreshold(2) & FixedTime & Valid)];
ByFit.Total = [sum(ByFit.Time),sum(ByFit.Distance)]; % Total Time & Distance cells in each trial type
ByFit.Ratio = [ByFit.Distance(1)/ByFit.Time(1),ByFit.Time(2)/ByFit.Distance(2)];
ByFit.MajIndex = [(ByFit.Distance(1)-ByFit.Time(1))/(ByFit.Distance(1)+ByFit.Time(1)),(ByFit.Distance(2)-ByFit.Time(2))/(ByFit.Time(2)+ByFit.Distance(2))];

disp(sum((ByFit.TimeCells.*ByFit.DistanceCells==1)));
disp(sum(ByFit.TimeCells)+sum(ByFit.DistanceCells));

% CellType metrics is based on the variance from the linear fit 
CellType = (R(:,20)-R(:,24))./(R(:,20)+R(:,24));
ByCellType.TimeCells= (CellType>CellTypeThreshold & Valid);
ByCellType.DistanceCells= (CellType<=-CellTypeThreshold & Valid);
ByCellType.Time=    [sum(CellType>CellTypeThreshold  & FixedDistance & Valid),sum(CellType>CellTypeThreshold  & FixedTime & Valid)];
ByCellType.Distance=[sum(CellType<=-CellTypeThreshold & FixedDistance & Valid),sum(CellType<=-CellTypeThreshold & FixedTime & Valid)];
ByCellType.Total = [sum(ByCellType.Time),sum(ByCellType.Distance)];
ByCellType.Ratio = [ByCellType.Distance(1)/ByCellType.Time(1),ByCellType.Time(2)/ByCellType.Distance(2)];
ByCellType.MajIndex = [(ByCellType.Distance(1)-ByCellType.Time(1))/(ByCellType.Distance(1)+ByCellType.Time(1)),(ByCellType.Distance(2)-ByCellType.Time(2))/(ByCellType.Distance(1)+ByCellType.Time(1))];

if ShowBar % Show the bar plot
    disp('PValue');
    disp(ByPValue);
    disp('Fit');
    disp(ByFit);
    disp('CellType');
    disp(ByCellType);
    
    figure(1);
    subplot(2,2,1);
    bFixedDistance = bar([ByCellType.Distance(1),ByCellType.Time(1);ByFit.Distance(1),ByFit.Time(1);ByPValue.Distance(1),ByPValue.Time(1)]);
    title('Fixed-Distance Sessions');
    xlabel(['CellType      ','           Fit              ','   P-Value']);
    xtips1 = bFixedDistance(1).XEndPoints;
    ytips1 = bFixedDistance(1).YEndPoints;
    labels1 = string(bFixedDistance(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    xtips2 = bFixedDistance(2).XEndPoints;
    ytips2 = bFixedDistance(2).YEndPoints;
    labels2 = string(bFixedDistance(2).YData);
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    
    subplot(2,2,2);
    bFixedDistance = bar([ByCellType.Distance(2),ByCellType.Time(2);ByFit.Distance(2),ByFit.Time(2);ByPValue.Distance(2),ByPValue.Time(2)]);
    title('Fixed-Time Sessions');
    xlabel(['CellType      ','           Fit              ','   P-Value']);
    xtips1 = bFixedDistance(1).XEndPoints;
    ytips1 = bFixedDistance(1).YEndPoints;
    labels1 = string(bFixedDistance(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    xtips2 = bFixedDistance(2).XEndPoints;
    ytips2 = bFixedDistance(2).YEndPoints;
    labels2 = string(bFixedDistance(2).YData);
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    
    subplot(2,2,3);
    bFixedDistancePer = bar([ByCellType.Distance(1)/(ByCellType.Distance(1)+ByCellType.Time(1)),ByCellType.Time(1)/(ByCellType.Distance(1)+ByCellType.Time(1));...
        ByFit.Distance(1)/(ByFit.Distance(1)+ByFit.Time(1)),ByFit.Time(1)/(ByFit.Distance(1)+ByFit.Time(1));ByPValue.Distance(1)/(ByPValue.Distance(1)+ByPValue.Time(1)),ByPValue.Time(1)/(ByPValue.Distance(1)+ByPValue.Time(1))]);
    xlabel(['CellType      ','           Fit              ','   P-Value']);
    ylim([0 1]);
    yticks([0:0.1:1]);
    grid on;
    
    subplot(2,2,4);
    bFixedDistancePer = bar([ByCellType.Distance(2)/(ByCellType.Distance(2)+ByCellType.Time(2)),ByCellType.Time(2)/(ByCellType.Distance(2)+ByCellType.Time(2));...
        ByFit.Distance(2)/(ByFit.Distance(2)+ByFit.Time(2)),ByFit.Time(2)/(ByFit.Distance(2)+ByFit.Time(2));ByPValue.Distance(2)/(ByPValue.Distance(2)+ByPValue.Time(2)),ByPValue.Time(2)/(ByPValue.Distance(2)+ByPValue.Time(2))]);
    xlabel(['CellType      ','           Fit              ','   P-Value']);
    ylim([0 1]);
    yticks([0:0.1:1]);
    grid on;
    
    sgtitle(file);
end

end