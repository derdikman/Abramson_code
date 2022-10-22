
% Analyze all session and save in summary.may, to be further analyzed in R and Excell
% Author : Shai Abramson Jan/2021

function RunAllSessionsTimeVsDistance


    % Fixed distance trials (525-900cm)
    R1 = NeuronFiringTimeVsDistance('bk26-0323',0);
    R2 = NeuronFiringTimeVsDistance('bk26-0326',0);
    R3 = NeuronFiringTimeVsDistance('bk26-0401',0);
    R4 = NeuronFiringTimeVsDistance('bk33-0723',0);
    R5 = NeuronFiringTimeVsDistance('bk35-0827',0);
    R6 = NeuronFiringTimeVsDistance('bk41-0317',0);
    R7 = NeuronFiringTimeVsDistance('bk41-0406',0);
    R8 = NeuronFiringTimeVsDistance('bk49-0217',0);
    R9 = NeuronFiringTimeVsDistance('bk49-0222',0);
    
    % Fixed time trials (16 sec)
    R10 = NeuronFiringTimeVsDistance('bk35-0831',0);
    R11 = NeuronFiringTimeVsDistance('bk35-0902',0);
    R12 = NeuronFiringTimeVsDistance('bk41-0325',0);
    R13 = NeuronFiringTimeVsDistance('bk41-0331',0);
    R14 = NeuronFiringTimeVsDistance('bk45-0803',0);
    R15 = NeuronFiringTimeVsDistance('bk45-0812',0);
    R16 = NeuronFiringTimeVsDistance('bk45-0826',0);
    R17 = NeuronFiringTimeVsDistance('bk49-0210',0);
    R18 = NeuronFiringTimeVsDistance('bk49-0214',0);

    save('Summary.mat','R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R17','R18');
    
end