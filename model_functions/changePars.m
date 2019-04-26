%% changePars


function csett = changePars(sett, loopId, l1, l2)
csett = struct;

csett.iIu = sett.iIu;
csett.eIu = sett.eIu;
csett.noiseamp = sett.noiseamp;
csett.Ipulse_start = sett.Ipulse_start;
csett.Ipulse_amp = sett.Ipulse_amp;
csett.Iper_freq = sett.Iper_freq;
csett.Iper_amp = sett.Iper_amp;
csett.g_syn_ee_r = sett.g_syn_ee_r;
csett.g_syn_ei_r = sett.g_syn_ei_r;
csett.g_syn_ii_r = sett.g_syn_ii_r;
csett.g_syn_ie_r = sett.g_syn_ie_r;
csett.lateonset_start = sett.lateonset_start;

for l = 1:sum(loopId)
    if l == 1
        ltemp = l1;
        Id = find(loopId > 0, 1, 'first');
    elseif l == 2
        ltemp = l2;
        Id = find(loopId > 0, 1, 'last');
    end
    
    switch Id
        case 1 % i current
            loopi_reg = find(sett.loopIi(1:sett.Nregions) == 1);
            
            if length(loopi_reg) == 1
                csett.iIu(loopi_reg) = sett.iIu(loopi_reg) + (ltemp-1)*sett.Istep;
            else
                csett.iIu(loopi_reg(l)) = sett.iIu(loopi_reg(l)) + (ltemp-1)*sett.Istep;
            end
            
        case 2 % e current
            loope_reg = find(sett.loopIe(1:sett.Nregions) == 1);
            
            if length(loope_reg) == 1
                csett.eIu(loope_reg) = sett.eIu(loope_reg) + (ltemp-1)*sett.Istep;
            else
                csett.eIu(loope_reg(l)) = sett.eIu(loope_reg(l)) + (ltemp-1)*sett.Istep;
            end
            
        case 3 % special current paths
            loopspec_reg = find(sett.loopIe(1:sett.Nregions) == 5);
            
            csett.iIu(loopspec_reg(l)) = sett.iIspecial{loopspec_reg(l)}(ltemp);
            csett.eIu(loopspec_reg(l)) = sett.eIspecial{loopspec_reg(l)}(ltemp);
            
        case 4 % noise
            csett.noiseamp = sett.noiseamp + (ltemp-1) * sett.noiseamp_step;
            
        case 5 % pulse
            if sett.Iapploop == 1 || (sett.Iapploop == 3 && l == 1)
                csett.Ipulse_start = sett.Ipulse_start + (ltemp-1)*sett.timestep;
            elseif sett.Iapploop == 2 || (sett.Iapploop == 3 && l == 2)
                csett.Ipulse_amp = sett.Ipulse_amp + (ltemp-1)*sett.ampstep;
            end
            
        case 6 % period current
            if sett.Iperloop == 1 || (sett.Iperloop == 3 && l == 1)
                csett.Iper_freq = sett.Iper_freq + (ltemp-1)*sett.perstep;
            elseif sett.Iperloop == 2 || (sett.Iperloop == 3 && l == 2)
                csett.Iper_amp = sett.Iper_amp + (ltemp-1)*sett.perampstep;
            end
            
        case 7 % interareal connections
            csett.g_syn_ei_r(sett.g_syn_loop == 1 | sett.g_syn_loop == 2) = (sett.g_syn_ei_r(sett.g_syn_loop == 1 | sett.g_syn_loop == 2) + sett.g_syn_ifactor * (ltemp-1)*sett.g_syn_r_step);
            csett.g_syn_ee_r(sett.g_syn_loop == 1 | sett.g_syn_loop == 2) = (sett.g_syn_ee_r(sett.g_syn_loop == 1 | sett.g_syn_loop == 2) + (ltemp-1)*sett.g_syn_r_step);% / sett.g_syn_ifactor;
            csett.g_syn_ie_r(sett.g_syn_loop == -1 | sett.g_syn_loop == 2) = (sett.g_syn_ie_r(sett.g_syn_loop == -1 | sett.g_syn_loop == 2) + sett.g_syn_ifactor * (ltemp-1)*sett.g_syn_r_step);
            csett.g_syn_ii_r(sett.g_syn_loop == -1 | sett.g_syn_loop == 2) = (sett.g_syn_ii_r(sett.g_syn_loop == -1 | sett.g_syn_loop == 2) + (ltemp-1)*sett.g_syn_r_step);% / sett.g_syn_ifactor;

        case 8 % delayed onset area 2
            csett.lateonset_start = sett.lateonset_start + (ltemp-1)*sett.lateonset_steps;
    end
end
end
