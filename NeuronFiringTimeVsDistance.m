
% Analyze a session data (filename) Neuron (Neuron#, 0=all)
% Generate figures and data summary (R)
% Author : Shai Abramson Jan/2021
%
% R legend :
% --------- 
% Total runs
% Valid runs
% Valid Runs Left Turns
% Valid Runs Right Turns
% Fit Time
% Fit Distance
% RMS
% Mean Time
% Mean Distance
% Max Onsset Distance
% Sum of Absolute residuals,
% R^2 (Distance vs Velocity)
% Pvalue (Distance vs Velocity)
% CellType based on Distance vs Velocity)
% Error of time (Distance vs Velocity)
% Error of distance (Distance vs Velocity)
% Error ot time, times Speed (Distance vs Velocity)
% R^2 (Time vs 1/Velocity)
% Pvalue (Time vs 1/Velocity)
% CellType (Time vs 1/Velocity)
% Error of Time (Time vs 1/Velocity)
% Error of Distance (Time vs 1/Velocity)
% Error of distance, times Speed (Time vs 1/Velocity)
% Err Slope
% R^2 Error
% P Error, 
% PeakFiringRate Low speed
% PeakFiringRate Medium speed
% PeakFiringRate High speed

%_______________________________________________________________________

function [R] = NeuronFiringTimeVsDistance(filename,Neuron)
warning('off');

tmp=load(strcat(filename,'.mat')); % Load the data
fn = fieldnames(tmp);
FirstVar = fn{1};
s = tmp.(FirstVar);

disp('')

TimeResolution = 0.2; % sec  , the analysis binning resolution
TimeResolutionHR = 0.1; % sec , the figures binning resolution

DistanceResolution = TimeResolution*20; % cm
DistanceResolutionHR = TimeResolutionHR*20; % cm

VelocityGroups = 3; % number of velocity groups
Vres = 5; % Velocity resolution on the raster image (10 == 0.1 cm/sec resolution)

% *************************************************************************************
minValidRuns = 10; % Minimum valid runs required for analysis
minPeak = 1; % Hz, min tuning curve peak firing frequency
minRate = 1; % Hz, min firing rate required for the 'first burst' search option
minSilence = 1; % Sec, min required silent period before onset
minPeakLocation = 2; % Sec, Generate figures only if tuning curve peak exceed this value
ExtraTime = 5; % Additional time to analyze after the treadmill stops
ShowPath = 0; % 0 -do not show path, 1 -Show the path
ShowPks = 0; % 0 -do not show, 1 -show the peaks values on the tuning curves
TurnType = 'B'; % The turn type to analyze, R(ight), L(eft) or B(both)
FigureFormat = 3; % 1- Including heat map, 2-Analyze 3-Paper 4- none

FigureType = '.bmp'; % use .epsc for high quality
OnsetType = 2; % 1- First burst, 2- minSilence before peak bin 3-max 
Shuffle = 0; % 0- no shuffle, 1 - Randomly shift cyclicly each run
simTimeCellSigma = 2; % Simulated time Cell sigma (sec)
Simulate = 0; % 0-none 1-simulate as Time cell 2- simulate as Distance cell

Titles = false; % Titles for figures
TrialTitle = false; % Information about the trail and neuron
FntSize = 14; % Font size
NumFntSize = 12; % numbers Font size


% *************************************************************************************


fprintf(2,'TM+%d minValidRuns=%d Format=%d OnsetType=%d minSilence=%d Simulate=%d \n',ExtraTime,minValidRuns,FigureFormat,OnsetType,minSilence,Simulate);

MxTime = 40; % Max time span
MxDistance = 1200;% Max distance span

MaxTime = round(MxTime/TimeResolution);  
MaxDistance = round(MxDistance/DistanceResolution);  

MaxTimeHR = round(MxTime/TimeResolutionHR);  
MaxDistanceHR = round(MxDistance/DistanceResolutionHR);  

NumberOfNeurons = size(s.neurons,1);
FiringNeuronCnt = 0; % Counter of Firing neurons
RejectedCells = 0; % Counter of rejected cells

if Neuron==0 % If Neuron parametr is 0, run through all
    nStart=1;
    nStop=NumberOfNeurons;
else % Run a specific neuron
    nStart=Neuron;
    nStop=Neuron;
end

% Iniitiate vectors
TRuns(1:nStop)=0;
VRuns(1:nStop)=0;
ValidRunsLeftTurns(1:nStop)=0;
ValidRunsRightTurns(1:nStop)=0;
FitT(1:nStop)=0;
FitD(1:nStop)=0;
RMS(1:nStop)=0;
MTime(1:nStop)=0;
MaxOnset(1:nStop)=0;
MDist(1:nStop)=0;
MaxOnsetDistance(1:nStop)=0;
NSumAbsResid(1:nStop)=0;
R2S_V(1:nStop)=0;
PS_V(1:nStop)=0;
CellS_V(1:nStop)=0;
Et_S_V(1:nStop)=0;
Es_S_V(1:nStop)=0;
Es_S_VxSpeed(1:nStop)=0;
R2T_RV(1:nStop)=0;
PT_RV(1:nStop)=0;
CellT_RV(1:nStop)=0;
Et_T_RV(1:nStop)=0;
Es_T_RV(1:nStop)=0;
Et_T_RVxSpeed(1:nStop)=0;
ErrSlope(1:nStop)=0;
R2Err(1:nStop)=0;
PErr(1:nStop)=0;
PeakFiringRateL(1:nStop)=0;
PeakFiringRateM(1:nStop)=0;
PeakFiringRateH(1:nStop)=0;

for n=nStart:nStop % Run through all neurons
    unit = s.neurons{n};
    ts = unit.timestamps ; % neurons firing timestamps
    Laptype = [s.zerothlap s.laptypes.']; % Forced direction after the treadmill (R or L)
    if TurnType=='B' % Both Right & Left turns
        tm.start = s.treadmill(:,1); % Treadmill start times
        tm.stop  = s.treadmill(:,2); % Treadmill stop times
        tm.Speed = s.treadmill(:,3); % Run speed
    else % Right or Left turns
        tm.start = s.treadmill(Laptype==TurnType,1); % Treadmill start times
        tm.stop  = s.treadmill(Laptype==TurnType,2); % Treadmill stop times
        tm.Speed = s.treadmill(Laptype==TurnType,3); % Run speed
    end
    
    RunPeriod =tm.stop-tm.start;
    if max(RunPeriod)-min(RunPeriod)>1 % If periods are the same, then a Time session, otherwise a Distance session
        SessionType=1; % Distanced session
    else
        SessionType=2; % Timed session
    end
       
    pT = s.path(:,1); % Path timestamp
    pX = s.path(:,2); % Path X coordinate
    pY = s.path(:,3);  % Path Y coordinate;
    pV = s.vel(:); % Path velocity
    
    if FigureFormat~=4 % If figures are reqested
        figure(1);
        clf;
        if TrialTitle
            sgtitle(strcat('Trial : ',filename,' Neuron# : ',int2str(n)),'Color','blue','FontSize',16);
        end
    end
    
    ptWithinTimeSlot = (pT>=0) & (pT<=1000000);  % During all runs
    % Interpolate position
    xs = interp1(pT(ptWithinTimeSlot),pX(ptWithinTimeSlot),ts,'nearest','extrap');
    ys = interp1(pT(ptWithinTimeSlot),pY(ptWithinTimeSlot),ts,'nearest','extrap');
    
    if FigureFormat == 1
        Map = zeros(round(1000/DistanceResolution), round(1000/DistanceResolution));  % Firing rate map
        for k=1:size(pX) % Path
            Mxs(k) = round(1+pX(k)/DistanceResolution);
            Mys(k) = round(1+pY(k)/DistanceResolution);
            Map(Mxs(k),Mys(k)) = 20;
        end
        for k=1:size(xs) % Heat map
            Mxs(k) = round(1+xs(k)/DistanceResolution);
            Mys(k) = round(1+ys(k)/DistanceResolution);
            Map(Mxs(k),Mys(k)) = Map(Mxs(k),Mys(k)) + 10;
        end
    end
    
    if FigureFormat~=4        
        % Draw the path of all trials
        figure(1);
        if ShowPath == 1
            ax = subplot(4,3,3); % Draw all trials path
            scatter(pX,pY,'.','black');
            title("\fontsize{"+num2str(FntSize)+"}Firing pos while and post treadmill");
            xlabel('\fontsize{'+num2str(FntSize)+'}Cm');
            ylabel('\fontsize{'+num2str(FntSize)+'}cm');
            hold on;
            
            ptWithinTimeSlot = (pT>=0) & (pT<=1000000);  % During all runs
            xs = interp1(pT(ptWithinTimeSlot),pX(ptWithinTimeSlot),ts,'nearest','extrap');
            ys = interp1(pT(ptWithinTimeSlot),pY(ptWithinTimeSlot),ts,'nearest','extrap');
            scatter(xs,ys,'.','red')
            axis([0 1000 0 1000]);
            
            ax = subplot(4,3,6); % Draw firing heat map
            image(flip(Map.',1));
            colormap(ax,'hot')
        end
    end
    
    Runs=size(tm.start,1);
    
    tmFiringTime(1:MaxTime)=0;
    tmFiringDistance(1:MaxDistance)=0;
    FiringRunCnt = 0;
    
    minTimedRun = round(min(tm.stop-tm.start)/TimeResolution+1);
    ShortestRunTime = min(tm.stop-tm.start);
    LongestRunTime = max(tm.stop-tm.start);
    minDistancedRun = round(min((tm.stop-tm.start).*tm.Speed)/DistanceResolution+1);
    ShortestRunDistance = min((tm.stop-tm.start).*tm.Speed);
    LongestRunDistance = max((tm.stop-tm.start).*tm.Speed);
    if Simulate==1 %Time cell simulation
        FieldCenter = normrnd(minTimedRun*0.5,5); % Normal distribution around the center of the minimum time trial
        if FieldCenter<1
            FieldCenter = 1;
        end
    end 
    if Simulate==2 % Distance cell simulation
        FieldCenter = normrnd(600,200); % Normal distribution around 600mm
        if FieldCenter<50
            FieldCenter = 50;
        end
    end
    for i=1:Runs  %Setup a vector with the firing rates on the treadmill, for each session
        start = tm.start(i);
        stop = tm.stop(i);
        AnalysisPeriod = stop-start+ExtraTime;
        tmWithinTimeSlot = (ts>=start) & (ts<=stop+ExtraTime);  
        tmFiringTime = (ts(tmWithinTimeSlot)-start); % Vector of firing times on the treadmill
        
        
        MaxCellRate = 5*max(accumarray(round(tmFiringTime/TimeResolution+1),1)); % Find the max rate for the simulation
        if MaxCellRate > 0
            if Simulate==1 % Time cell simulation
                tmFiringTime = SimTimeCell(FieldCenter,simTimeCellSigma,MaxCellRate,60).'; %Simulate time cell with FieldCenter, simTimeCellSigma sigma and MaxCellrate for 60sec total trial
            end
            if Simulate==2 % Distance cell simulation
                tmFiringTime = SimDistanceCell(FieldCenter/tm.Speed(i),3,MaxCellRate,60).'; %Simulate distance cell with FieldCenter, 3sec sigma and MaxCellrate for 60sec total trial
            end
        end
        
        tmFiringTime = tmFiringTime(tmFiringTime<AnalysisPeriod);
        if size(tmFiringTime,1) == 0
            tmFiringTime(1)=0.1;
        end
              
        if Shuffle == 1
            Shift = random('Uniform',0,AnalysisPeriod); % Time shift in case of Shuffle, randomed uniformly from the range 0:AnalysisPeriod
            tmFiringTime(i) = mod(tmFiringTime + Shift, AnalysisPeriod);
        end
        
        tmAllFiringTime{i} = tmFiringTime;
        tmFiringDistance = tmFiringTime*tm.Speed(i); % vector of firing distances on the treadmill     

        tmFiringRatePerTime=accumarray(round(tmFiringTime/TimeResolution+1),1); % bining by time
        tmFiringRatePerTimeHR=accumarray(round(tmFiringTime/TimeResolutionHR+1),1); % bining by time, high resolution
        
        tmFiringRatePerDistance=accumarray(round(tmFiringDistance/DistanceResolution+1),1); % bining by distance
        tmFiringRatePerDistanceHR=accumarray(round(tmFiringDistance/DistanceResolutionHR+1),1); % bining by distance, high resolution
        
        tmRunSizeTime=size(tmFiringRatePerTime,1);
        tmRunSizeDistance=size(tmFiringRatePerDistance,1);
        tmRunSizeTimeHR=size(tmFiringRatePerTimeHR,1);
        tmRunSizeDistanceHR=size(tmFiringRatePerDistanceHR,1);
        
        % Save all vectors for this neuron
        FiringRunCnt = FiringRunCnt + 1;
        AllFiringTime(FiringRunCnt,:) = [tmFiringRatePerTime;zeros(MaxTime-tmRunSizeTime,1)]; % Collect all firing time vectors
        AllFiringDistance(FiringRunCnt,:) = [tmFiringRatePerDistance;zeros(MaxDistance-tmRunSizeDistance,1)]; % Collect all firing distance vectors
        
        AllFiringTimeHR(FiringRunCnt,:) = [tmFiringRatePerTimeHR;zeros(MaxTimeHR-tmRunSizeTimeHR,1)]; % Collect all firing time vectors
        AllFiringDistanceHR(FiringRunCnt,:) = [tmFiringRatePerDistanceHR;zeros(MaxDistanceHR-tmRunSizeDistanceHR,1)]; % Collect all firing distance vectors
               
        % Color the positions of the firings on the TM and after the TM
        if FigureFormat~=4
            figure(1);
            ax = subplot(4,3,3);
        end
        if ShowPath == 1
            ptWithinTimeSlot = (pT>=start) & (pT<=stop); % Draw (in Green) the firing positions while on the treadmill
            xs = interp1(pT(ptWithinTimeSlot),pX(ptWithinTimeSlot),ts,'nearest','extrap');
            ys = interp1(pT(ptWithinTimeSlot),pY(ptWithinTimeSlot),ts,'nearest','extrap');
            scatter(xs,ys,'.','green')
            axis([0 1000 0 1000]);
            
            ptWithinTimeSlot = (pT>=stop) & (pT<=stop+ExtraTime); % Draw (in Yellow) the firing positions from treadmill stop until ExtraTime
            xs = interp1(pT(ptWithinTimeSlot),pX(ptWithinTimeSlot),ts,'nearest','extrap');
            ys = interp1(pT(ptWithinTimeSlot),pY(ptWithinTimeSlot),ts,'nearest','extrap');
            scatter(xs,ys,'.','yellow')
            axis([0 1000 0 1000]);
        end
    end
      
    hold off;
    
    [~,SpeedOrder] = sort(tm.Speed); % Sort the tm speeds vector
    
    SortedAllFiringTime = AllFiringTime(SpeedOrder,:); % sort the firing vectors according to the speeds vector
    SortedAllFiringDistance = AllFiringDistance(SpeedOrder,:);
    
    SortedAllFiringTimeHR = AllFiringTimeHR(SpeedOrder,:); % sort the firing vectors according to the speeds vector
    SortedAllFiringDistanceHR = AllFiringDistanceHR(SpeedOrder,:);
    
    Onset = SortedAllFiringTime>=minRate; % Map of all bursts exceeding minRate threshold
    FirstOnset = zeros(size(Onset));
         
    % Calculate tuning curve per speed group
    for u=1:VelocityGroups
        NormalizationFactor = VelocityGroups/Runs/TimeResolution; % The averaging factor (number of runs in every velocity group) times the bin width (Time Resolution)
        FiringTimePerSpeedGroups1(u,:) = smoothdata(sum(SortedAllFiringTime(round(1+(u-1)/VelocityGroups*Runs):round(u*Runs/VelocityGroups),:))*NormalizationFactor);
        FiringTimePerSpeedGroups(u,:) = smoothdata(FiringTimePerSpeedGroups1(u,:));
        FiringDistancePerSpeedGroups1(u,:) = smoothdata(sum(SortedAllFiringDistance(round(1+(u-1)/VelocityGroups*Runs):round(u*Runs/VelocityGroups),:))*NormalizationFactor);
        FiringDistancePerSpeedGroups(u,:) = smoothdata(FiringDistancePerSpeedGroups1(u,:));
    end
    
    % Choose to find the peaks per the session type (this is particulary
    % aimed to find the post-TM peak and since it is usually a Place cell
    % at the corridor the burst will best align per the session type
    if SessionType == 1 % Fixed-Distance
        pp1 = FiringDistancePerSpeedGroups(1,:);
        pp2 = FiringDistancePerSpeedGroups(2,:);
        pp3 = FiringDistancePerSpeedGroups(3,:);
    else % Fixed-time
        pp1 = FiringTimePerSpeedGroups(1,:);
        pp2 = FiringTimePerSpeedGroups(2,:);
        pp3 = FiringTimePerSpeedGroups(3,:);
    end
      
    pp1(1)=0; % Set the first point to 0 so that findpeaks can detect a peak at the start
    % Find bin peak with minimum minPeak , min width of 3 bins (width is
    % at 3db points) and a seperation of 10 bins 
    [pks1,plocs1] = findpeaks(pp1,'MinPeakHeight',minPeak,'MinPeakWidth',3,'MinPeakDistance',2,'NPeaks',2);
    Npks1 = numel(pks1); % Number of peaks
    

    pp2(1)=0; % Set the first point to 0 so that findpeaks can detect a peak at the start
    [pks2,plocs2] = findpeaks(pp2,'MinPeakHeight',minPeak,'MinPeakWidth',3,'MinPeakDistance',2,'NPeaks',2);
    Npks2 = numel(pks2); % Number of peaks
    

    pp3(1)=0; % Set the first point to 0 so that findpeaks can detect a peak at the start
    [pks3,plocs3] = findpeaks(pp3,'MinPeakHeight',minPeak,'MinPeakWidth',3,'MinPeakDistance',2,'NPeaks',2);
    Npks3 = numel(pks3); % Number of peaks
    
    Npks = max([Npks1,Npks2,Npks3]); % Take the max # of peaks from the 3 velocity groups
    plocs(1)=9999;
    if Npks ~= 0
        if Npks1~=0
            plocs(1) = plocs1(1);
        end
        if Npks2~=0
            plocs(1) = min(plocs2(1),plocs(1));
        end
        if Npks3~=0
            plocs(1) = min(plocs3(1),plocs(1));
        end
    end
    plocs(2)=9999;
    if Npks==2
        if Npks1==2
            plocs(2) = plocs1(2);
        end
        if Npks2==2
            plocs(2) = min(plocs2(2),plocs(2));
        end
        if Npks3==2
            plocs(2) = min(plocs3(2),plocs(2));
        end
    end

    % Find the peak firing rate in each velocity group
    [PeakFiringRate(1), PeakFiringLocation(1)]=max(FiringTimePerSpeedGroups(1,:));
    [PeakFiringRate(2), PeakFiringLocation(2)]=max(FiringTimePerSpeedGroups(2,:));
    [PeakFiringRate(3), PeakFiringLocation(3)]=max(FiringTimePerSpeedGroups(3,:));
    
    minPeakLocationOfAllGroups = min(PeakFiringLocation);
    
    if (max(FiringTimePerSpeedGroups(:))>minPeak || max(FiringDistancePerSpeedGroups(:))>minPeak) && minPeakLocationOfAllGroups>minPeakLocation
        % Calculate peaks, their positions, center of mass and width of tuning curve
        for u=1:VelocityGroups
            [PeakT(n,u), PosT(n,u)] = max(FiringTimePerSpeedGroups(u,:)); % Peaks and positions of the tuning curves
            [PeakD(n,u), PosD(n,u)] = max(FiringDistancePerSpeedGroups(u,:));
            
            % Calculate center of mass
            CM = 0;
            for j=1:size(FiringTimePerSpeedGroups(u,:),2)
                CM = CM + j*FiringTimePerSpeedGroups(u,j);
            end
            CMassT(n,u) = CM/sum(FiringTimePerSpeedGroups(u,:));
            CM = 0;
            for j=1:size(FiringDistancePerSpeedGroups(u,:),2)
                CM = CM + j*FiringDistancePerSpeedGroups(u,j);
            end
            CMassD(n,u) = CM/sum(FiringDistancePerSpeedGroups(u,:));
            
            % Calculate width
            W = find(FiringTimePerSpeedGroups(u,:)>(PeakT(u)/2),1,'Last')-find(FiringTimePerSpeedGroups(u,:)>(PeakT(u)/2),1,'First')+1; % width of the tuning curve
            if size(W,2)==0
                W = 0;
            end
            WidthT(n,u) = W;
            W = find(FiringDistancePerSpeedGroups(u,:)>(PeakD(u)/2),1,'Last')-find(FiringDistancePerSpeedGroups(u,:)>(PeakD(u)/2),1,'First')+1; % width of the tuning curve
            if size(W,2)==0
                W = 0;
            end
            WidthD(n,u) = W;
        end
        
        % Find correlation between the tuning curves
        for u=1:VelocityGroups-1 
            [PCorT(n,u), CT] = max(conv(FiringTimePerSpeedGroups(u,:),fliplr(FiringTimePerSpeedGroups(u+1,:))));
            [PCorD(n,u), CD] = max(conv(FiringDistancePerSpeedGroups(u,:),fliplr(FiringDistancePerSpeedGroups(u+1,:))));
            CorT(n,u) = CT - size(FiringTimePerSpeedGroups,2);
            CorD(n,u) = CD - size(FiringDistancePerSpeedGroups,2);
        end
        
        % Prepare the search range for finding the onset position       
        RejectCell=false;
        if Simulate == 0
            MaxTimeToAnalyzePostTM = ExtraTime; % By default analyze to the end of TM + ExtraTime
        else
            MaxTimeToAnalyzePostTM = 0; % On cells simulations analyze to the end of TM
        end
        if Npks == 2 % In case of 2 peaks, reject cell if 2nd peak is post TM
            if (SessionType == 1 && plocs(2)*DistanceResolution>round(ShortestRunDistance+40)) || (SessionType == 2 && plocs(2)*TimeResolution>round(ShortestRunTime+1)) % If peak is post shortest TM+1sec/40cm then reject cell
                RejectCell = true;
                RejectedCells = RejectedCells + 1;
            end
        else if Npks == 1 % In case of single peak post TM, reject the cell
                if (SessionType == 1 && plocs(1)*DistanceResolution>round(ShortestRunDistance+40)) || (SessionType == 2 && plocs(1)*TimeResolution>round(ShortestRunTime+1)) % If peak is post shortest TM+1sec/40cm then reject cell
                    RejectCell = true;
                    RejectedCells = RejectedCells + 1;
                end
            end
        end
        
        % Find the Onset Position  
        OnSetPos(1:Runs) = 0;
        for i=1:Runs
            if MaxTimeToAnalyzePostTM>0
                MaxTimeToAnalyze = tm.stop(i)-tm.start(i)+ MaxTimeToAnalyzePostTM;
            else
                MaxTimeToAnalyze = tm.stop(i)-tm.start(i);
            end
            AnalyzedTime(i) = MaxTimeToAnalyze;
            AllFiringTimeOnTM = AllFiringTime(i,1:round(MaxTimeToAnalyze/TimeResolution)); % a vector containing only the firing time on the treadmill
            if OnsetType == 3 % Find the max (not really onset, actually  max)
                [~,MaxPos] = max(AllFiringTimeOnTM.'); % Find the peak firing on the treadmill
                if MaxPos > 1
                    OnSetPos(i)=MaxPos;
                    OnSetTime(i)= MaxPos*TimeResolution;
                    FirstOnset(i,OnSetPos(i)) = 2;
                end
            else
                if OnsetType == 2 % Find first null before peak
                    [~,MaxPos] = max(AllFiringTimeOnTM.'); % Find the peak firing on the treadmill
                    if MaxPos > 1
                        OnSetPos(i)=MaxPos;
                        Spikes = tmAllFiringTime{i}(tmAllFiringTime{i}<(MaxPos+1)*TimeResolution); % Collect all spikes timestamps before the peak bin
                        OnSetTime(i)= Spikes(end);
                        for k = size(Spikes,1):-1:2
                            Gap = Spikes(k)-Spikes(k-1);
                            if Gap <minSilence % if gap between spikes <minSilence then set the Onset to the earlier spike
                                OnSetTime(i)= Spikes(k-1);
                            else
                                break; % If gap >minSilence then stop searching
                            end
                        end
                        OnSetPos(i) = round(OnSetTime(i)/TimeResolution+1); % Onset Bin
                        FirstOnset(i,OnSetPos(i)) = 2;
                    end
                else % OnsetType=1. First burst
                    OPos = find(AllFiringTimeOnTM>=minRate);
                    if size(OPos,2) >0
                        OnSetPos(i) = OPos(1); % The first burst exceeding minRate
                        FirstOnset(i,OnSetPos(i)) = 2;
                    else
                        OnSetPos(i) = 0;
                    end
                end
            end
        end
        if size(OnSetPos,2)~=Runs
            fprintf('error n=%d \n',n);
            ValidRuns = 0;
        else
            % Collect only the points with the valid Onset for the linear regression
            ValidRuns = 0;
            RejectedRuns = 0;
            clear Rx; % Valid runs index
            clear Ry; % Valid runs Onset time
            clear Rv; % Valid runs speeds
            clear RTurns; % Valid runs Turn type (1=Left 0=Right)
            if (Npks>0) && ~RejectCell % Analyze Onsets only if a peak was detected in the tuning curve
                for k = 1:Runs
                    if (OnSetPos(k)~=0) 
                        ValidRuns = ValidRuns+1;
                        Rx(ValidRuns)=k; % Store only the valid runs positions
                        Rv(ValidRuns)=tm.Speed(k); % Store the valid runs speeds
                        Rt(ValidRuns)=OnSetTime(k); % Store the valid runs onset time (sec)
                        Rd(ValidRuns)=Rt(ValidRuns)*Rv(ValidRuns); % Store the valid runs onset distance (cm)
                        RTurns(ValidRuns)= (Laptype(k)=='L'); % Store the valid runs turn type
                    end
                end
            end
        end
        
        TRuns(n) = Runs;    % Total Runs
        VRuns(n) = ValidRuns/Runs; % Ratio of Valid Runs
        if ValidRuns<minValidRuns % Ignore trials with less than minValidRuns
            fprintf('%s N=%d ValidRuns=%d (< %d) Rejected=%d \n',filename,n,ValidRuns,minValidRuns, RejectCell);
        else      

            R = [1:Runs];
            
            V = tm.Speed(Rx).'; % Velocities of the valid runs
            [SortedV,Vorder] = sort(V); % Sorted velocities
            Vround = round(SortedV*Vres).'; % Rounded velocities, Vres resolution
            
            Time = OnSetTime(Rx); % Onset time
            Distance = V.*Time; % Onset distance
            ConstTime = mean(Time)*[0:MxTime]; % Distance in case of constant time
            ConstTime1 = mean(Time).*ones(1,(MxTime+1)); % Constant time
            ConstDist = mean(Distance).*ones(1,(MxTime+1)); % Constant distance
            minSpeed = min(tm.Speed);
            maxSpeed = max(tm.Speed);
            if SessionType==1 % Fixed Distance
                minTimeValue = min(tm.stop-tm.start); % min time
                minTime = minTimeValue.*ones(size(minSpeed:maxSpeed)); % min time vector
                minDist = minTimeValue.*[minSpeed:maxSpeed]; % distance of the min time
            else % Fixed Time
                minDistValue = min((tm.stop(SpeedOrder)-tm.start(SpeedOrder)).*tm.Speed(SpeedOrder)); % min distance
                minDist = minDistValue.*ones(size(minSpeed:maxSpeed)); % min Dist vector
                minTime = minDistValue./[minSpeed:maxSpeed]; % min distance
            end
            
            ConstDist1 = mean(Distance)./V; % Time in case of constant distance
            clear V_fit;
            
            % D(v)= v*t+D0
            % A perfect distance cell at distance D : D(V)=D
            % A perfect time cell at time T : D(V)=V*T
            O = ~isoutlier(Distance,'gesd'); % Non outlier elements of Distance
            [fitDistTime,S] = polyfit(Time(O),Distance(O) ,1); % Linear regression between the Distance at the OnSet and its time
            [RegDistTime,] = polyval(fitDistTime,Time,S); % Estimation of distance based on time
            [RegDistTime1,] = polyval(fitDistTime,[0:MxTime],S); % Estimation of the distance based on time, entire range
            mdlS_T = fitlm(Time,Distance); 
            
            [fitTimeRVel,S] = polyfit(1./V(O),Time(O),1); % Linear regression between the OnSet time and 1/velocity
            [RegTimeRVel,] = polyval(fitTimeRVel,1./V,1); % Estimation of time based on 1/velocity
            [RegTimeRVel1,] = polyval(fitTimeRVel,1./[0:MxTime],1); % Estimation of time based on 1/velocity, entire range
            [RegTimeRVel2,] = polyval(fitTimeRVel,1./tm.Speed(SpeedOrder),1); % Estimation of time based on 1/velocity, entire range
            
            [r2 rmse] = rsquare(Time,RegTimeRVel,false);
            [h,p] = ttest(RegTimeRVel,Time);
            [cR,cP]=corrcoef(Time,1./V);
            cS = size(cR);
            if cS(1)==2
                cR=cR(2,1);
            end
            mdlRV_T = fitlm(Time,1./V);
            mdlT_RV = fitlm(1./V,Time);
            eval(strcat(filename(1:4),filename(6:end),'{',int2str(n),'}.TVR=mdlT_RV;'));          
            
            m = mdlT_RV.Coefficients.Estimate(2); % regression coefficient
            m0 = mdlT_RV.Coefficients.Estimate(1); % regression constant
            
            EsT_RV = sum((m./V-Time+m0).^2);
            EtT_RV = sum((Time-mean(Time)).^2);
            EtT_RVxSpeed = sum(((Time-mean(Time)).*V).^2);

            CellTypeT_RV = (EtT_RV-EsT_RV)/(EtT_RV+EsT_RV); % -1 for Time Cell, +1 for Place Cell
                    
            [fitDistVel,S] = polyfit(V(O),Distance(O),1); % Linear regression between the OnSet distance and the velocity
            [RegDistVel,] = polyval(fitDistVel,V,1); % Estimation of distance based on velocity
            [RegDistVel1,] = polyval(fitDistVel,[0:MxTime],1); % Estimation of distance based on velocity, entire range
            [DVr2 DVrmse] = rsquare(Distance,RegDistVel,false);
            [DVh,DVp] = ttest(RegDistVel,Distance);
            [CDVr,CDVp]=corrcoef(Distance,V);
            CDVs = size(CDVr);
            if CDVs(1)==2
                CDVr=CDVr(2,1);
                CDVp=CDVp(2,1);
            end
            mdlS_V = fitlm(V,Distance);
            eval(strcat(filename(1:4),filename(6:end),'{',int2str(n),'}.SV=mdlS_V;'));

            m = mdlS_V.Coefficients.Estimate(2); % regression coefficient
            m0 = mdlS_V.Coefficients.Estimate(1); % regression constant
            EtS_V = sum((m*V-Distance+m0).^2); % Residuals
            EsS_V = sum((Distance-mean(Distance)).^2); % Residuals
            EsS_VxSpeed = sum(((Distance-mean(Distance))./V).^2); % Residuals

            CellTypeS_V = (EtS_V-EsS_V)/(EtS_V+EsS_V); % -1 for Time Cell, +1 for Place Cell          
            
            % Sort the residual errors by runs order and make a linear regression
            TimeResid = abs(RegTimeRVel-OnSetPos(Rx));
            DistanceResid = abs(RegDistVel-OnSetPos(Rx).*V);
            [fitTimeE,~] = polyfit(Rx,TimeResid,1); % Linear regression between the Run (corresponding to velocity) and the onset time
            [TimeE_fit,] = polyval(fitTimeE,Rx,S); %Regression error evaluation
            [fitDistanceE,~] = polyfit(Rx,DistanceResid,1); % Linear regression between the Run (corresponding to velocity) and the onset time
            [DistanceE_fit,] = polyval(fitDistanceE,Rx,S); %Regression error evaluation
            
            [ResidR2 Residrmse] = rsquare(TimeResid,TimeE_fit,false);
            [h1,p1] = ttest(TimeE_fit,TimeResid,'Alpha',0.01);
            [cR1,cP1]=corrcoef(TimeE_fit,Rx);
            cS1 = size(cR1);
            if cS1(1)==2
                cR1=cR1(2,1);
            end
            mdlTimeErr = fitlm(TimeResid,Rx);
            mdlDistanceErr = fitlm(DistanceResid,Rx);
            
            FitT(n) = fitDistVel(1); %  Regression coefficient of Distance vs Velocity
            FitD(n) = fitTimeRVel(1); %  Regression coefficient of Time vs 1/Velocity
            RMS(n) = S.normr/sqrt(ValidRuns); % RMS error of Onset distance estimation
            MTime(n) = mean(Time); % Average Onset time
            [M,MPos] = max(OnSetPos); % Maximum OnSet 
            MaxOnset(n) = M*TimeResolution; % Maximum OnSet time
            MaxOnsetDistance(n) = MaxOnset(n)*tm.Speed(MPos); % Maximum OnSet Distance
            MDist(n) = mean(Distance); % Average Onset distance
            NSumAbsResid(n) = sum(DistanceResid) / ValidRuns; % Average Distance residual
            ErrSlope(n) = fitDistanceE(1); % Regression slope of the residuals sorted by Runs
            ValidRunsLeftTurns(n) = sum(RTurns); % Valid runs with Left turns after TM
            ValidRunsRightTurns(n) = ValidRuns-sum(RTurns); % Valid runs with Right turns after TM
            
            fprintf('%s N=%d Runs= %d Valid=%d Left=%d Right=%d \n',filename,n, TRuns(n),round(VRuns(n)*TRuns(n)),ValidRunsLeftTurns(n),ValidRunsRightTurns(n));
            
            R2S_V(n) = round(mdlS_V.Rsquared.Ordinary,2); % R2
            PS_V(n) = round(mdlS_V.Coefficients.pValue(2),2); % P value
            Es_S_V(n) = EsS_V;
            Et_S_V(n) = EtS_V;            
            Es_S_VxSpeed(n) = EsS_VxSpeed;
            CellS_V(n) = round(CellTypeS_V,2); % 
            PeakFiringRateL(n)=PeakFiringRate(1);
            PeakFiringRateM(n)=PeakFiringRate(2);
            PeakFiringRateH(n)=PeakFiringRate(3);
            
            R2T_RV(n) = round(mdlT_RV.Rsquared.Ordinary,2); % R2
            PT_RV(n) = round(mdlT_RV.Coefficients.pValue(2),2); % P value
            CellT_RV(n) = round(CellTypeT_RV,2); % Cell Type based on Time vs 1/Velocity
            Et_T_RV(n) = EtT_RV;
            Es_T_RV(n) = EsT_RV;
            Et_T_RVxSpeed(n) = EtT_RVxSpeed;
           
            R2Err(n) = round(mdlDistanceErr.Rsquared.Ordinary,2); % R2
            PErr(n) = round(mdlDistanceErr.Coefficients.pValue(2),2); % P value
                       
            FiringTimeImage = (SortedAllFiringTimeHR>0)*1;
            FiringDistanceImage = (SortedAllFiringDistanceHR>0)*1;
            
            for j=1:Runs  % Mark the different velocity groups
                if j>=Runs*2/3
                    FiringTimeImage(j,:) = FiringTimeImage(j,:)*3;
                    FiringDistanceImage(j,:) = FiringDistanceImage(j,:)*3;
                else
                    if j>=Runs/3
                        FiringTimeImage(j,:) = FiringTimeImage(j,:)*2;
                        FiringDistanceImage(j,:) = FiringDistanceImage(j,:)*2;
                    end
                end
            end
            
            for j=1:ValidRuns  % Mark the Onset point on the valid runs
                ixt = round(Rt(j)/TimeResolutionHR+1);
                ixd = round(Rd(j)/DistanceResolutionHR+1);
                ix = find(SpeedOrder==Rx(j)); % Position of the respective valid vector in the sorted velocity array
                FiringTimeImage(ix,ixt) = 4;
                FiringDistanceImage(ix,ixd) = 4;
            end
            FiringTimeImage = FiringTimeImage+1; % Backgound >0
            FiringDistanceImage = FiringDistanceImage+1; % Backgound >0
                       
            switch FigureFormat
                case 1
                    ax = subplot(4,3,1);
                case 2
                    ax = subplot(4,2,1);
                case 3
                    ax = subplot(4,1,1);
            end
            
            if FigureFormat<=3
                image([0 MxTime],[],FiringTimeImage);
                cmap = [ 0.9 0.9 0.9; 0 0 0.7 ; 0.5 0 0; 0 0.5 0; 0 0.5 0];
                colormap(ax, cmap);
                axis([0 40 0 Runs]);
                ax = gca;
                ax.FontSize = NumFntSize;
                ax.YDir = 'reverse';
                box(ax,'off')
%                 title("\fontsize{"+num2str(FntSize)+"}Firing rate per time, sorted by Velocity");
                if Titles
                    ylabel("\fontsize{"+num2str(FntSize)+"}Sorted Run#"+newline);
                    xlabel("\fontsize{"+num2str(FntSize)+"}Sec"+newline);
                end
                ax.FontSize = NumFntSize;
                ax.XTick = [0 40];
                ax.YTick = [0 Runs];
                hold on;
                plot((tm.stop(SpeedOrder)-tm.start(SpeedOrder)),(1:Runs),'color',[0.2 0.2 0.2],'LineStyle','--','LineWidth',1); % The treadmill end time
                plot((tm.stop(SpeedOrder)-tm.start(SpeedOrder)+ExtraTime),(1:Runs),'color',[0 0 0],'LineStyle',':','LineWidth',1); % The treadmill end time + ExtraTime
                hold off;
            end
            
            switch FigureFormat
                case 1
                    ax = subplot(4,3,2);
                case 2
                    ax = subplot(4,2,2);
                case 3
                    ax = subplot(4,1,1);
            end
                       
            if FigureFormat<3
                image([0 MxDistance],[],FiringDistanceImage);
                cmap = [ 1 1 1; 0 0 0.9 ; 0.9 0 0; 1 0.9 0; 0 0 0];
                colormap(ax, cmap);
                ax = gca;
                ax.FontSize = NumFntSize;
                ax.YDir = 'reverse';
                box(ax,'off')
                if Titles
                    ylabel("\fontsize{"+num2str(FntSize)+"}Sorted Run#"+newline);
                    xlabel("\fontsize{"+num2str(FntSize)+"}Cm"+newline);
                end
                hold on;
                plot((tm.stop(SpeedOrder)-tm.start(SpeedOrder)).*tm.Speed(SpeedOrder),(1:Runs),'green'); % The treadmill end time
                hold off;
            end
     
            switch FigureFormat
                case 1
                    ax = subplot(4,3,4);
                case 2
                    ax = subplot(4,2,3);
                case 3
                    if SessionType==2
                        ax = subplot(4,1,2);
                    else
                        ax = subplot(4,1,2);
                    end
            end
            
            if FigureFormat<=3
                plot(FiringTimePerSpeedGroups.','LineWidth',3);
                set(gca,'ColorOrder',[0 0 0.5; 0.5 0 0; 0 0.5 0]);
                axis([0 80 0 round(1+max(FiringTimePerSpeedGroups(:)))]);
                ax = gca;
                ax.FontSize = NumFntSize;
                box(ax,'off');
                ax.XTick = [0 80];
                ax.XTickLabels = [0 40];
                ax.YTick = [0 round(1+max(FiringTimePerSpeedGroups(:)))];
                if Titles
                    ylabel("\fontsize{"+num2str(FntSize)+"}Firing rate (Hz)"+newline);
                    xlabel("\fontsize{"+num2str(FntSize)+"}Sec"+newline);
                end
                if SessionType ==2 && FigureFormat~=2 && ShowPks ==1
                    text(plocs1(:)+10,pks1(:),strcat(num2str(pks1(:),2),' Hz'),'color','blue');
                    text(plocs2(:)+10,pks2(:),strcat(num2str(pks2(:),2),' Hz'),'color','red');
                    text(plocs3(:)+10,pks3(:),strcat(num2str(pks3(:),2),' Hz'),'color','yellow');
                end
            end
                        
            switch FigureFormat
                case 1
                    ax = subplot(4,3,5);
                case 2
                    ax = subplot(4,2,4);
                case 3
                    if SessionType==2
                        ax = subplot(4,1,3);
                    else
                        ax = subplot(4,1,3);
                    end
            end
            
            if FigureFormat<=3
                plot(FiringDistancePerSpeedGroups.','LineWidth',3);
                set(gca,'ColorOrder',[0 0 0.5; 0.5 0 0; 0 0.5 0]);
                axis([0 120 0 round(1+max(FiringDistancePerSpeedGroups(:)))]);
                ax.YTick = [0 round(1+max(FiringDistancePerSpeedGroups(:)))];
                set(gca,'XTick',[0 120]);
                set(gca,'XTickLabels',[0 1200],'FontSize',FntSize);
                ax = gca;
                ax.FontSize = NumFntSize;
                box(ax,'off')
                if Titles
                    xlabel("\fontsize{"+num2str(FntSize)+"}Cm"+newline);
                    ylabel("\fontsize{"+num2str(FntSize)+"}Firing rate (Hz)"+newline);
                end
                if SessionType ==1 && FigureFormat~=2 && ShowPks ==1
                    text(plocs1(:)+10,pks1(:),strcat(num2str(pks1(:),2),' Hz'),'color','blue');
                    text(plocs2(:)+10,pks2(:),strcat(num2str(pks2(:),2),' Hz'),'color','red');
                    text(plocs3(:)+10,pks3(:),strcat(num2str(pks3(:),2),' Hz'),'color','yellow');
                end
            end
            
            switch FigureFormat
                case 1
                    ax = subplot(4,3,7);
                case 2
                    ax = subplot(4,2,5);
                case 3
                    ax = subplot(4,1,4);
            end
                       
            FirstOnsetByVelocity = zeros(max(Vround)-min(Vround)+1,size(FirstOnset,2));
            
            for k = 1:ValidRuns
                indx = round(V(k)*Vres)-min(Vround)+1; % index by the velocity of the k Run
                FirstOnsetByVelocity(indx,:)= FirstOnsetByVelocity(indx,:) + FirstOnset(Rx(k),:);
            end
            
            if FigureFormat<=3
                cla;
                image([0 MxTime],[min(Vround/Vres) max(Vround/Vres)],2*FirstOnsetByVelocity(:,:)); % Raster map of Onset bin sorted by 1/velocities order
                cmap = [ 1 1 1; 0.9 0.9 0.9 ; 1 0 0];
                colormap(ax, cmap);
                axis([0 MxTime floor(min(Vround/Vres)) ceil(max(Vround/Vres)+1)]);
                ax = gca;
                ax.FontSize = NumFntSize;
                ax.YDir = 'reverse';
                box(ax,'off')
                if Titles
                    xlabel("\fontsize{"+num2str(FntSize)+"}Sec"+newline);
                    ylabel("\fontsize{"+num2str(FntSize)+"}Velocity (Cm/Sec)"+newline);
                end
                ax.XTick = [0 MxTime];
                ax.XTickLabels = [0 MxTime];
                ax.YTick = [floor(min(Vround/Vres)) ceil(max(Vround/Vres)+1)];
                hold on;
                plot(200*Vres./Vround,Vround/Vres,'blue'); % A curve representing the time to readh 200cm distance for each run
                plot(600*Vres./Vround,Vround/Vres,'blue'); % A curve representing the time to readh 600cm distance for each run
                plot(1000*Vres./Vround,Vround/Vres,'blue'); % A curve representing the time to readh 1000cm distance for each run
                plot(1400*Vres./Vround,Vround/Vres,'blue'); % A curve representing the time to readh 1400cm distance for each run
                plot(1800*Vres./Vround,Vround/Vres,'blue'); % A curve representing the time to readh 1800cm distance for each run
                plot((tm.stop(SpeedOrder)-tm.start(SpeedOrder)),tm.Speed(SpeedOrder),'color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1); % The treadmill end time
                plot(AnalyzedTime(SpeedOrder),tm.Speed(SpeedOrder),'color',[0.3 0.3 0.3],'LineStyle',':'); % The treadmill end time + ExtraTime
                                
                [RegTimeRVel3,] = polyval(fitTimeRVel,Vres./Vround,1); % Estimation of time based on 1/velocity, entire range
                plot(RegTimeRVel3,Vround/Vres,'black'); % Estimated Onset based on regression
                hold off;
            end
            
            switch FigureFormat
                case 1
                    ax = subplot(4,3,8);
                case 2
                    ax = subplot(4,2,6);
            end
                        
            if FigureFormat<3
                image([0 MxDistance],[],0*FirstOnset); %
                cmap = [ 1 1 1; 0.9 0.9 0.9 ; 1 0 0];
                colormap(ax, cmap);
                ax = gca;
                ax.FontSize = NumFntSize;
                box(ax,'off')
                xlabel("\fontsize{"+num2str(FntSize)+"}Cm"+newline);
                ylabel("\fontsize{"+num2str(FntSize)+"}Run#"+newline);
                hold on;
                scatter(DistanceResid,Rx,'.','red');
                plot(DistanceE_fit,Rx,'black','LineWidth',1);
                hold off;
            end
            
            switch FigureFormat
                case 1
                    ax = subplot(4,3,10);
                case 2
                    ax = subplot(4,2,7);
            end

            if FigureFormat<3               
                scatter(Rv,Distance,'x','black');
                ax = gca;
                ax.FontSize = NumFntSize;
                box(ax,'off')
                if Titles
                    xlabel("\fontsize{"+num2str(FntSize)+"}Velocity (Cm/sec)");
                    ylabel("\fontsize{"+num2str(FntSize)+"}Distance (Cm)");
                end
                hold on;
                %             scatter(Rv(~O),Distance(~O),'x','red'); % show outliers in red
                axis([25 50 -inf inf]);
                plot((0:MxTime),ConstTime,'blue','LineWidth',1); %
                plot((0:MxTime),ConstDist,'red','LineWidth',1);
                plot((0:MxTime),RegDistVel1,'black','LineWidth',2);
                plot(tm.Speed(SpeedOrder),(tm.stop(SpeedOrder)-tm.start(SpeedOrder)).*tm.Speed(SpeedOrder),'green'); % The treadmill end time
                hold off;
            end
            
            switch FigureFormat
                case 1
                    ax = subplot(4,3,11);
                case 2
                    ax = subplot(4,2,8);
            end
                       
            if FigureFormat<3
                hold on;
                scatter(1./V,Time,'x','black');
                ax = gca;
                ax.FontSize = NumFntSize;
                box(ax,'off')
                xlabel("\fontsize{"+num2str(FntSize)+"}1/Velocity (Sec/Cm)");
                ylabel("\fontsize{"+num2str(FntSize)+"}Time (Sec)");
                hold on;
                axis([0.02 0.04 -inf inf]);
                plot(1./[0:MxTime],ConstTime1,'blue','LineWidth',1); %
                plot(1./V,ConstDist1,'red','LineWidth',1);
                plot(1./V,RegTimeRVel,'black','LineWidth',2);
                plot(1./tm.Speed(SpeedOrder),(tm.stop(SpeedOrder)-tm.start(SpeedOrder)),'green'); % The treadmill end time
                hold off;
            end
            
            if FigureFormat~=4 % Save only figures with significant pValue 
                PlaceCell = mdlT_RV.Coefficients.pValue(2)<0.05;
                TimeCell = mdlS_V.Coefficients.pValue(2)<0.05;
                Both = TimeCell & PlaceCell;
                if Both
                    saveas(figure(1),strcat('MBoth-',filename,'_N',int2str(n),'#2',FigureType));
                else if TimeCell  % Time Cell
                        saveas(figure(1),strcat('MTime-',filename,'_N',int2str(n),'#2',FigureType));
                    else if PlaceCell  % Place Cell
                            saveas(figure(1),strcat('MPlace=',filename,'_N',int2str(n),'#2',FigureType));
                        else
                            saveas(figure(1),strcat('M',filename,'_N',int2str(n),'#2',FigureType));
                        end
                    end
                end
            end
            
        end
    end
end

fprintf('Rejected Cells = %d \n',RejectedCells);
R = [TRuns.',VRuns.',ValidRunsLeftTurns.',ValidRunsRightTurns.',FitT.', FitD.', RMS.',MTime.',MaxOnset.',MDist.',MaxOnsetDistance.',NSumAbsResid.',R2S_V.',PS_V.',CellS_V.',Et_S_V.',Es_S_V.',Es_S_VxSpeed.',R2T_RV.',PT_RV.',CellT_RV.',Et_T_RV.',Es_T_RV.',Et_T_RVxSpeed.',ErrSlope.',R2Err.',PErr.',PeakFiringRateL.',PeakFiringRateM.',PeakFiringRateH.'];
save(strcat('mdls_',filename),strcat(filename(1:4),filename(6:end)));


function [r2 rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).
if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end
% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end
% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);
if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end
rmse = sqrt(mean((y(:) - f(:)).^2));
