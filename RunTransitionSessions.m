
% Analyze all session and save in summary.may, to be further analyzed in R and Excell
% Author : Shai Abramson Jan/2021

function RunTransitionSession

    % Fixed distance trials (525-900cm)

    R7 = NeuronFiringTimeVsDistance('bk41-0406',0); % Transition from fixed distance to fixed time
    R8 = NeuronFiringTimeVsDistance('bk49-0217',0); % Transition from fixed distance to fixed time
    
    % Fixed time trials (16 sec)
    R10 = NeuronFiringTimeVsDistance('bk35-0831',0); % Transition from fixed time to fixed distance
    R12 = NeuronFiringTimeVsDistance('bk41-0325',0); % Transition from fixed time to fixed distance


    save('SummaryTransitionSession.mat','R7','R8','R10','R12');
    
end