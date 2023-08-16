
% The script plots the averaged max firing rate of each run vs the
% velocity, in Fixed-Time and Fixed-Distance experiments

close all
clear global

Dir='G:\Dropbox (Personal)\Brain\Project';
Animals = {'BK49','BK41','BK35','BK45','BK33','BK26'};
NumOfAnimals = 6; % Use 1 only for a short validity test, otherwise 5
ThetaCalc = "median"; % 'median' or 'average'
BinSize = 1; %Sec since we are looking for the max firing in every run, we should bin at relatively large bins
VBin=2; % Velocity bin size Cm/sec

cnt=0;
for AnimalNum=1:NumOfAnimals
    trial=[];
    
    switch AnimalNum
        case 1
            Animal = '\KrausData\Recordings\BK49';
            trial{1} = 'bk49-0210';
            trial{2} = 'bk49-0214';
            trial{3} = 'bk49-0217';
            trial{4} = 'bk49-0222';
            trial_type = [1 1 2 2]; % 1=fixed-time, 2=fixed-distance
            
        case 2
            Animal = '\KrausData\Recordings\BK41';
            trial{1} = 'bk41-0325';
            trial{2} = 'bk41-0331';
            trial{3} = 'bk41-0317';
            trial{4} = 'bk41-0406';
            trial_type = [1 1 2 2]; % 1=fixed-time, 2=fixed-distance
            
        case 3
            Animal = '\KrausData\Recordings\BK35';
            trial{1} = 'bk35-0831';
            trial{2} = 'bk35-0902';
            trial{3} = 'bk35-0827';
            trial_type = [1 1 2]; % 1=fixed-time, 2=fixed-distance
            
        case 4
            Animal = '\KrausData\Recordings\BK45';
            trial{1} = 'bk45-0803';
            trial{2} = 'bk45-0812';
            trial{3} = 'bk45-0826';
            trial_type = [1 1 1]; % 1=fixed-time, 2=fixed-distance
            
        case 5
            Animal = '\KrausData\Recordings\BK33';
            trial{1} = 'bk33-0723';
            trial_type = [2]; % 1=fixed-time, 2=fixed-distance
            
        case 6
            Animal = '\KrausData\Recordings\BK26';
            trial{1} = 'bk26-0323';
            trial{2} = 'bk26-0326';
            trial{3} = 'bk26-0401';
            trial_type = [2 2 2]; % 1=fixed-time, 2=fixed-distance
    end
    
    
    
    %---------------------------------------------------------------------
    
    for Tr=1:length(trial)
        d=[];
        teta=[];
        load([Dir '\' trial{Tr}]);
        VarName = eval([trial{Tr}(1:4) '_' trial{Tr}(end-3:end)]);
        TM= VarName.alltreadmill;
        disp([AnimalNum, Tr,length(VarName.neurons)]);
        for Neuron=1:length(VarName.neurons)
            for j=1:length(TM) % Go through all runs within the trial
                ts= VarName.neurons{Neuron}.timestamps;
                ts=ts(ts>=TM(j,1) & ts<=TM(j,2));
                ts = ts-min(ts);
                % Define the bin edges
                if ~isempty(ts)
                    binEdges = 0:BinSize:TM(j,4);
                    [counts, binCenters] = histcounts(ts, binEdges);
                    Vel = TM(j,3);
                    mx=max(counts/BinSize);
                    if mx<=15 % Count only neurons with max firing rate<15hz (other are inter-neurons or noise)
                        cnt=cnt+1;
                        maxF(cnt,1)=mx; % max Firing rate
                        maxF(cnt,2)=Vel; % Run Speed
                        maxF(cnt,3)=trial_type(Tr); % Trial Type
                    end
                end
            end
        end
    end
end

f1=figure(1);
clf;

FTime = maxF(:,3)==1; % Fixed-time data
FDist = maxF(:,3)==2; % Fixed-distance data

meanFiring_at_V=[];
stdErrFiring_at_V=[];
vcnt=0;
T= 25:VBin:(50-VBin);
for V=25:VBin:(50-VBin) 
    vcnt=vcnt+1;
    FTimeV = FTime & maxF(:,2)>=V & (maxF(:,2)<=V+VBin);
    meanFiring_at_V(vcnt)=mean(maxF(FTimeV,1));
    stdErrFiring_at_V(vcnt)=std(maxF(FTimeV,1))/sqrt(numel(maxF(FTimeV,1)));
end
b1=bar(T, meanFiring_at_V,0.35,'FaceColor','b');
hold on;
errorbar(T, meanFiring_at_V, stdErrFiring_at_V, 'k.');

meanFiring_at_V=[];
stdErrFiring_at_V=[];
vcnt=0;
for V=25:VBin:(50-VBin) 
    vcnt=vcnt+1;
    FDistV = FDist & maxF(:,2)>=V & (maxF(:,2)<=V+VBin);
    meanFiring_at_V(vcnt)=mean(maxF(FDistV,1));
    stdErrFiring_at_V(vcnt)=std(maxF(FDistV,1))/sqrt(numel(maxF(FDistV,1)));
end
b2=bar(T+0.7, meanFiring_at_V,0.35,'FaceColor','r');
hold on;
errorbar(T+0.7, meanFiring_at_V, stdErrFiring_at_V, 'k.');


xlabel('Velocity cm/sec');
ylabel('Average Max firing rate (Hz)');
xlim([24,49]);

title('Average Max firing rate vs velocity at Fixed-Time and Fixed-Distance experiments');
legend([b1,b2],{'Fixed Time','Fixed-Distance'});
