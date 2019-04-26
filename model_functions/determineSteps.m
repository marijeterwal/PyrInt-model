%% determine l1steps and l2steps

function [out] = determineSteps(sett, Id, nr)

switch Id
    case 1 % i current   
        loopi_reg = find(sett.loopIi(1:sett.Nregions) == 1);
        out = round((sett.iIu_max(loopi_reg(1:nr)) - sett.iIu(loopi_reg(1:nr))) / sett.Istep) + 1;   
        
    case 2 % e current
        loope_reg = find(sett.loopIe(1:sett.Nregions) == 1);
        out = round((sett.eIu_max(loope_reg(1:nr)) - sett.eIu(loope_reg(1:nr))) / sett.Istep) + 1; 
        
    case 3 % special current paths
        loopspec_reg = find(sett.loopIe(1:sett.Nregions) == 5);
        for n = 1:nr
            out(n) = length(sett.eIspecial{loopspec_reg(n)}); 
        end
        if sum(ismember(loopspec_reg,find(sett.loopIi(1:sett.Nregions) == 5))) < length(loopspec_reg)
            warning('Special loop settings for i and e are not identical, please check!')
        end
        
    case 4 % noise
        out = round((settnoiseamp_max - sett.noiseamp) / sett.noiseamp_step) + 1;
        
    case 5 % pulse
        if sett.Iapploop == 1 || sett.Iapploop == 3
            out(1) = round((sett.Ipulse_end - sett.Ipulse_start) / sett.timestep) +1;
        end
        if sett.Iapploop == 2 || sett.Iapploop == 3
            out(nr) = round((sett.Ipulse_amp_max - sett.Ipulse_amp) / sett.ampstep) +1;
        end
        
    case 6 % period current  
        if sett.Iperloop == 1 || sett.Iperloop == 3
        out(1) = round((sett.Iper_freq_max - sett.Iper_freq) / sett.perstep) +1;
        end
        if sett.Iperloop == 2 || sett.Iperloop == 3
            out(nr) = round((sett.Iper_amp_max - sett.Iper_amp) / sett.perampstep) +1;
        end
        
    case 7 % interareal connections
            if sum(sett.g_syn_loop(:))<0
                out = round((sett.g_syn_r_max(1) - sett.g_syn_ee_r(sett.g_syn_loop == -1)) / sett.g_syn_r_step) + 1;
            else
                out = round((sett.g_syn_r_max(1) - sett.g_syn_ee_r(sett.g_syn_loop == 1 | sett.g_syn_loop == 2)) / sett.g_syn_r_step) + 1;
            end
%         
    case 8 % delayed onset
        out = round((sett.lateonset_end - sett.lateonset_start) / sett.lateonset_steps) +1;
end
    
