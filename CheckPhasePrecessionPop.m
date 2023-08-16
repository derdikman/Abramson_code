
% Check the Time and Distance cells for a given Phase Precession slope PVal or CorrCoef
% PPSummary provides the mean values of the slope and R2 for a given pValue (<=)
% Then count the cells, for each type and experiment type, complying to the given mean slope criteria

function CheckPhasePrecessionPop

ConditionType=2; % 1 for pValue, 2 for CorrCoefficient
dSlopeTh = 0;
tSlopeTh = 0;
PValueTh = 0.05;
CorrCoefTh = -0.2;
PPSummary = SummarizePhasePrecession('AnalyzedData_Peak+5sec_200msBins');

TimeCells = [PPSummary.ByCellType]==1;
DistanceCells = [PPSummary.ByCellType]==2;

for i=1:length(PPSummary)
    if isempty(PPSummary(i).tSlope)
        tSlope(i) = 0;
    else
        tSlope(i) = PPSummary(i).tSlope;
    end
    if isempty(PPSummary(i).dSlope)
        dSlope(i) = 0;
    else
        dSlope(i) = PPSummary(i).dSlope;
    end
    if isempty(PPSummary(i).tPValue)
        tPValue(i) = 1;
    else
        tPValue(i) = PPSummary(i).tPValue;
    end
    if isempty(PPSummary(i).dPValue)
        dPValue(i) = 1;
    else
        dPValue(i) = PPSummary(i).dPValue;
    end
    if isempty(PPSummary(i).tCorr)
        tCorr(i) = 0;
    else
        tCorr(i) = PPSummary(i).tCorr;
    end
    if isempty(PPSummary(i).dCorr)
        dCorr(i) = 0;
    else
        dCorr(i) = PPSummary(i).dCorr;
    end
end

if ConditionType==1
    timeSlope = tSlope<tSlopeTh & tPValue<PValueTh; % Check the Slope criteria based on the time pp slope
    distSlope = dSlope<dSlopeTh & dPValue<PValueTh; % Check the Slope criteria based on the distance pp slope
else
    timeSlope = tCorr<CorrCoefTh; % Check the Correlation Coefficient criteria
    distSlope = dCorr<CorrCoefTh; % Check the Correlation Coefficient criteria
end
% Fixed Time
T_FT=[sum(TimeCells([PPSummary.Trial]>9 & timeSlope)),sum(TimeCells([PPSummary.Trial]>9 & distSlope)),sum(TimeCells([PPSummary.Trial]>9))];
D_FT=[sum(DistanceCells([PPSummary.Trial]>9 & timeSlope)),sum(DistanceCells([PPSummary.Trial]>9 & distSlope)),sum(DistanceCells([PPSummary.Trial]>9))];

% Fixed Distance
T_FD=[sum(TimeCells([PPSummary.Trial]<10 & timeSlope)),sum(TimeCells([PPSummary.Trial]<10 & distSlope)),sum(TimeCells([PPSummary.Trial]<10))];
D_FD=[sum(DistanceCells([PPSummary.Trial]<10 & timeSlope)),sum(DistanceCells([PPSummary.Trial]<10 & distSlope)),sum(DistanceCells([PPSummary.Trial]<10))];

Totals = ([D_FD(3) D_FT(3); T_FD(3) ,T_FT(3)]);
tdata = round([D_FD(1) D_FT(1); T_FD(1) ,T_FT(1)]./[D_FD(3) D_FT(3); T_FD(3) ,T_FT(3)],2);
ddata = round([D_FD(2) D_FT(2); T_FD(2) ,T_FT(2)]./[D_FD(3) D_FT(3); T_FD(3) ,T_FT(3)],2);
leg = {'FDist', 'FTime'; 'Dist', 'Time'};  % Example legend
T0 = table(Totals(1,:)', Totals(2,:)', 'VariableNames', leg(1,:), 'RowNames', leg(2,:));
T1 = table(tdata(1,:)' , tdata(2,:)' , 'VariableNames', leg(1,:), 'RowNames', leg(2,:));
T2 = table(ddata(1,:)' , ddata(2,:)' , 'VariableNames', leg(1,:), 'RowNames', leg(2,:));



figure(1);
clf;
subplot(2,1,1);
plot([PPSummary(TimeCells & [PPSummary.Trial]>9).tSlope],[PPSummary(TimeCells & [PPSummary.Trial]>9).tR2],'.b');
hold on;
plot([PPSummary(DistanceCells & [PPSummary.Trial]>9).tSlope],[PPSummary(DistanceCells & [PPSummary.Trial]>9).tR2],'.r');
plot([PPSummary(TimeCells & [PPSummary.Trial]<10).tSlope],[PPSummary(TimeCells & [PPSummary.Trial]<10).tR2],'ob');
plot([PPSummary(DistanceCells & [PPSummary.Trial]<10).tSlope],[PPSummary(DistanceCells & [PPSummary.Trial]<10).tR2],'or');
xlim([-300 0]);
xlabel('Phase Precession rate deg/sec');
ylabel('R2');
legend('Time Cells & FixedTime','Distance Cells & FixedTime','Time Cells & FixedDistance','Distance Cells & FixedDistance');
xticks([-300:20:0]);
set (gca,'Xdir','reverse');

grid on;

subplot(2,1,2);
plot([PPSummary(TimeCells & [PPSummary.Trial]>9).dSlope],[PPSummary(TimeCells & [PPSummary.Trial]>9).dR2],'.b');
hold on;
plot([PPSummary(DistanceCells & [PPSummary.Trial]>9).dSlope],[PPSummary(DistanceCells & [PPSummary.Trial]>9).dR2],'.r');
plot([PPSummary(TimeCells & [PPSummary.Trial]<10).dSlope],[PPSummary(TimeCells & [PPSummary.Trial]<10).dR2],'ob');
plot([PPSummary(DistanceCells & [PPSummary.Trial]<10).dSlope],[PPSummary(DistanceCells & [PPSummary.Trial]<10).dR2],'or');
xlim([-5 0]);
xticks([-5:1:0]);
xlabel('Phase Precession rate deg/cm');
ylabel('R2');
set (gca,'Xdir','reverse');
grid on;
% legend('Time Cells & FixedTime','Distance Cells & FixedTime','Time Cells & FixedDistance','Distance Cells & FixedDistance');

if ConditionType==1
    fprintf(2,['\nPValue<',num2str(PValueTh),' dSlope<',num2str(dSlopeTh),' tSlope<',num2str(tSlopeTh),'\n\n']);
else
    fprintf(2,['\nCorrCoef<',num2str(CorrCoefTh),' dSlope<',num2str(dSlopeTh),' tSlope<',num2str(tSlopeTh),'\n\n']);
end
disp('Total cells');
disp(T0)
disp('% of cells with pp rate<0 deg/sec');
disp(T1)
disp('% of cells with pp rate<0 deg/cm');
disp(T2)

end