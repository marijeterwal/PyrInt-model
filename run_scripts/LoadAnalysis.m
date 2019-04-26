% load analysis

if runnumber < 650%  %? 652
    load([sett.savelocdata,'analysis', '.mat']);
    if reg > 1;
        load([sett.savelocdata,'analysis2', '.mat']);
%         load([sett.savelocdata,'analysis3', '.mat']);
%         load([sett.savelocdata,'analysis4', '.mat']);
%         load([sett.savelocdata,'analysis5', '.mat']);
    end
    
else
    load([sett.savelocdata,'analysis', '.mat']);
    if reg > 1;  load([sett.savelocdata,'analysis2', '.mat']); end
    if loopId(5) > 0;  load([sett.savelocdata,'analysis_pulses', '.mat']); end
    if analyzeInfo; load([sett.savelocdata,'analysis_info', '.mat']); end
    if analyzeDirection; load([sett.savelocdata,'analysis_direction', '.mat']); end
end