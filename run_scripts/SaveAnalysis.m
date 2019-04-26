% saveAnalysis

if savedata
    save([sett.savelocdata,'analysis', '.mat'], varNames_standard1{:}); 
    if reg > 1;  save([sett.savelocdata,'analysis2', '.mat'], varNames_standard2{:}); end
    if loopId(5) > 0;  save([sett.savelocdata,'analysis_pulses', '.mat'], varNames_pulses{:}); end
    if analyzeInfo; save([sett.savelocdata,'analysis_info', '.mat'], varNames_info{:}); end
end