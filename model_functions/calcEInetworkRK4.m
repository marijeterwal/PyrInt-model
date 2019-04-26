%% INFO
% all matrices in this code use the format #neurons  x #time steps, since
% accessing a single column is faster than accessing a single row

function [spikesi, spikese] = calcEInetworkRK4(sett, para, csett)

Ni = sett.Ni;
Ne = sett.Ne;
dt = sett.dt;
reg = sett.Nregions;
tnumb = sett.Ttot/dt +1;

%% Parameters - to speed up computation in for loop
iphi = para.iphi;

%Reverse potentials
iE_l = para.iE_l;
iE_na = para.iE_na;
iE_k = para.iE_k;

eE_l = para.eE_l;
eE_na = para.eE_na;
eE_k = para.eE_k;

E_ampa = para.E_ampa;
E_gaba = para.E_gaba;

%Conductances
ig_l = para.ig_l;
ig_na = para.ig_na;
ig_k = para.ig_k;
eg_l = para.eg_l;
eg_na = para.eg_na;
eg_nap = para.eg_nap;
eg_kdr = para.eg_kdr;
eg_ka = para.eg_ka;

%Synapses
itheta_syn = para.itheta_syn;
etheta_syn = para.etheta_syn;
ialpha_syn = para.ialpha_syn;
ibeta_syn = para.ibeta_syn;
ealpha_syn = para.ealpha_syn;
ebeta_syn = para.ebeta_syn;

%Others
eC_m = para.eC_m;
iC_m = para.iC_m;

thresholdi      = 0;
thresholde      = -20;

%% Settings and currents
seed            = sett.seed;
fixed_initcond  = sett.fixed_initcond;
delay           = sett.delay;
delay_reg_e     = sett.delay_reg_e;
delay_reg_i     = sett.delay_reg_i;

% Looped settings:
names = fieldnames(csett);
for n = 1:length(names)
    sett.(names{n}) = csett.(names{n});
end

[conee_l, conei_l, conii_l, conie_l, conee_r, conei_r, conii_r, conie_r] = connectionmatrix(sett, 1);
[Ii, Ie, Ii_dyn, Ie_dyn, ipulsevct, epulsevct] = currents(sett, para);

% late onset
Iisave = Ii;
Iesave = Ie;  

%% Pulses
% pulse - these mimick synaptic inputs and are therefore voltage-dependent.
% The current value is therefore calculated in the RK loop
if (sum(sett.Iipulse) + sum(sett.Iepulse) >= 1) 
    pulseind = 1;
else
    pulseind = 0;
end
Iipulse = zeros(Ni*reg, length(sett.Ipulse_start));
Iepulse = zeros(Ne*reg, length(sett.Ipulse_start));
Ii_pulse = zeros(Ni*reg, 1);
Ie_pulse = zeros(Ne*reg, 1);

pulse_start = sett.Ipulse_start;

%% preallocate memory

lengths = delay_reg_i/dt+1;

Vi = zeros(Ni*reg, 2);
si = zeros(Ni*reg, lengths);

Ve = zeros(Ne*reg, 2);
se = zeros(Ne*reg, lengths);

KVi = zeros(Ni*reg, 4);
Kni = zeros(Ni*reg, 4);
Khi = zeros(Ni*reg, 4);
Ksi = zeros(Ni*reg, 4);

KVe = zeros(Ne*reg, 4);
Kbe = zeros(Ne*reg, 4);
Khe = zeros(Ne*reg, 4);
Kne = zeros(Ne*reg, 4);
Kse = zeros(Ne*reg, 4);

tempspikese = cell(tnumb,1);
tempspikesi = cell(tnumb,1);
spiketimese = cell(tnumb,1);
spiketimesi = cell(tnumb,1);


%% initial conditions
s = RandStream('mcg16807','Seed',seed);
RandStream.setGlobalStream(s);
if fixed_initcond == 2
    Vi(:,2) = -20*rand(Ni*reg,1)-50; %uniformly distributed between -50 and -70 mV
elseif fixed_initcond == 1
    Vi(:,2) = -65;
end

s = RandStream('mcg16807','Seed',seed+1);
RandStream.setGlobalStream(s);
if fixed_initcond == 2
    Ve(:,2) = -20*rand(Ne*reg,1)-50;
elseif fixed_initcond == 1
    Ve(:,2) = -75;
end

%interneurons
alphai_h0 = 0.07*exp(-(Vi(:,2)+58)/20); %steady state values for h, n and s
betai_h0 = 1./(exp(-0.1*(Vi(:,2)+28))+1);
alphai_n0 = -0.01*(Vi(:,2)+34)./(exp(-0.1*(Vi(:,2)+34))-1);
betai_n0 = 0.125* exp(-(Vi(:,2)+44)/80);
Fi0 = 1./(1 + exp(-(Vi(:,2)-itheta_syn)/2));
hi = alphai_h0./(alphai_h0 + betai_h0);
ni = alphai_n0./(alphai_n0 + betai_n0);
si(:,2) = ialpha_syn.* Fi0./(ialpha_syn*Fi0 + ibeta_syn);

%pyramidal neurons
be = 1 ./ (exp( (Ve(:,2) + 80) / 6) + 1); %steady state values for b, h, n, z and s
he = 1 ./ (exp( (Ve(:,2) + 53) / 7) + 1);
ne = 1 ./ (exp(-(Ve(:,2)+ 30) / 10) + 1);
Fe0 = 1./(1 + exp(-(Ve(:,2)-etheta_syn)/2));
se(:,2) = ealpha_syn.* Fe0./(ealpha_syn*Fe0 + ebeta_syn);

%% Seeding for noise
s = RandStream('mcg16807','Seed',seed+30);
RandStream.setGlobalStream(s);

%% Run RK4

Wt2 = [0 .5 .5 1];
rkIndex = [1 1 2 3];

dtiC = dt * (1/iC_m);
dtiphi = dt * iphi;
dteC = dt * (1/eC_m);
sdt = sqrt(dt);
siIlambda = sqrt(6*sett.iIlambda);
seIlambda = sqrt(6*sett.eIlambda);

for T = 2:tnumb
    
    Tsprev = mod(T-1, lengths) + 1;
    Tsnew = mod(T, lengths) + 1;
    Tvprev = mod(T-1, 2) + 1;
    Tvnew = mod(T, 2) + 1;
    
    %% Synaptic delays
    if T <= delay/dt
        si_delay = zeros(Ni*reg,1);
        se_delay = zeros(Ne*reg,1);
    else
        si_delay = si(:,mod(T -(delay/dt), lengths) +1);
        se_delay = se(:,mod(T -(delay/dt), lengths) +1);
    end
    
    if reg > 1
        if T <= delay_reg_e/dt
            se_delay_reg_e = zeros(Ne*reg,1);
            si_delay_reg_e = zeros(Ni*reg,1);
        else
            se_delay_reg_e = se(:,mod(T -(delay_reg_e/dt), lengths) +1);
            si_delay_reg_e = si(:,mod(T -(delay_reg_e/dt), lengths) +1);
        end

        if T <= delay_reg_i/dt
            se_delay_reg_i = zeros(Ne*reg,1);
            si_delay_reg_i = zeros(Ni*reg,1);
        else
            se_delay_reg_i = se(:,mod(T -(delay_reg_i/dt), lengths) +1);
            si_delay_reg_i = si(:,mod(T -(delay_reg_i/dt), lengths) +1);
        end
    end
    
    Sie = conie_l * si_delay;
    Sii = conii_l * si_delay;
    See = conee_l * se_delay;
    Sei = conei_l * se_delay;
    if reg > 1
        Sei_r = conei_r * se_delay_reg_i;
        See_r = conee_r * se_delay_reg_e;
        Sii_r = conii_r * si_delay_reg_i;
        Sie_r = conie_r * si_delay_reg_e;
    end
    
    %% late onset
    if sett.lateonset == true && T <= sett.lateonset_start/dt
            Ii(Ni+1:2*Ni) = zeros(Ni,1);
            Ie(Ne+1:2*Ne) = zeros(Ne,1);
    elseif sett.lateonset == true && T > sett.lateonset_start/dt
        Ii = Iisave;
        Ie = Iesave;
    end

    %% Runge-Kutta
    
    for rk = 1:4  %Fourth Order Runge-Kutta
        
        %definitions
        Vih = Vi(:,Tvprev) + KVi(:, rkIndex(rk))*Wt2(rk);
        Veh = Ve(:,Tvprev) + KVe(:, rkIndex(rk))*Wt2(rk);

        sih = si(:,Tsprev) + Ksi(:, rkIndex(rk))*Wt2(rk);
        seh = se(:,Tsprev) + Kse(:, rkIndex(rk))*Wt2(rk);
        
        hih = hi + Khi(:, rkIndex(rk))*Wt2(rk);
        nih = ni + Kni(:, rkIndex(rk))*Wt2(rk);
        beh = be + Kbe(:, rkIndex(rk))*Wt2(rk);
        heh = he + Khe(:, rkIndex(rk))*Wt2(rk);
        neh = ne + Kne(:, rkIndex(rk))*Wt2(rk);
        
        %interneurons
        alpha_m = (-0.1*Vih-3.5)./(exp(-0.1*Vih-3.5)-1); % divide by zero is possible
        beta_m = 4*exp(-0.0556*Vih-3.3333);
        alpha_h = 0.07*exp(-0.05*Vih-2.9);
        beta_h = 1./(exp(-0.1*Vih-2.8)+1);
        alpha_n = (-0.01*Vih-0.34)./(exp(-0.1*Vih-3.4)-1); % divide by zero is possible
        beta_n = 0.125* exp(-0.0125*Vih-0.55);
        
        alpha_m(abs(Vih + 35) < 0.001) = 1; % solving divide by zero problem
        alpha_n(abs(Vih + 34) < 0.001) = 1;
        
        mi_inf = alpha_m./(alpha_m + beta_m);
        
        % Pyramidal neurons
        me_inf = 1 ./ (exp(-0.1053*Veh-3.1579)+1);
        he_inf = 1 ./ (exp(0.1429*Veh+7.5714)+1);
        pe_inf = 1 ./ (exp(-0.20*Veh-8)+1);
        ne_inf = 1 ./ (exp(-0.10*Veh-3) +1);
        ae_inf = 1 ./ (exp(-0.05*Veh-2.5) +1);
        be_inf = 1 ./ (exp(0.1667*Veh+13.3333)+1);
        the = 0.37 + 2.78 ./ (exp(0.1667*Veh+6.75)+1);
        tne = 0.37 + 1.85 ./ (exp(0.0667*Veh+1.8)+1);
        
        %synapses
        Fi = 1./(1 + exp(-0.5*(Vih-itheta_syn)));
        Fe = 1./(1 + exp(-0.5*(Veh-etheta_syn)));
        
        %Currents
        Ii_l = ig_l * (Vih-iE_l);
        Ii_na = ig_na * mi_inf .* mi_inf .* mi_inf .* hih .* (Vih-iE_na);
        Ii_k = ig_k * nih .* nih .* nih .* nih .* (Vih-iE_k);
        
        Ie_l   = eg_l   * (Veh - eE_l);
        Ie_na  = eg_na  * me_inf .* me_inf .* me_inf  .* heh .* (Veh - eE_na);
        Ie_nap = eg_nap * pe_inf .* (Veh - eE_na);
        Ie_kdr = eg_kdr * neh .* neh .* neh .* neh .* (Veh - eE_k);
        Ie_ka  = eg_ka * ae_inf .* ae_inf .* ae_inf .* beh .* (Veh - eE_k);
        
        %Synaptic currents
        I_syn_ii = Sii .* (Vih - E_gaba);
        I_syn_ie = Sie .* (Veh - E_gaba);
        I_syn_ee = See .* (Veh - E_ampa);
        I_syn_ei = Sei .* (Vih - E_ampa);
        if reg > 1
            I_syn_ee_r = See_r .* (Veh - E_ampa);
            I_syn_ei_r = Sei_r .* (Vih - E_ampa);
            I_syn_ie_r = Sie_r .* (Veh - E_gaba);
            I_syn_ii_r = Sii_r .* (Vih - E_gaba);
        else
            I_syn_ee_r = zeros(Ne*reg,1);
            I_syn_ei_r = zeros(Ni*reg,1);
            I_syn_ie_r = zeros(Ne*reg,1);
            I_syn_ii_r = zeros(Ni*reg,1);
        end        

        %Pulses
        if pulseind == 1
        for k = 1:length(pulse_start)
            if T*dt >= pulse_start(k) && T*dt < (pulse_start(k) + 30)
                Iipulse(:,k) = ipulsevct(:,(T - round(pulse_start(k)/dt))+1) .* (Vih - E_ampa);
                Iepulse(:,k) = epulsevct(:,(T - round(pulse_start(k)/dt))+1) .* (Veh - E_ampa);
            else
                Iipulse(:,k) = zeros(Ni*reg, 1);
                Iepulse(:,k) = zeros(Ne*reg, 1);
            end
        end
        Ii_pulse = sum(Iipulse, 2);
        Ie_pulse = sum(Iepulse, 2);
        end

        % K-values
        KVi(:,rk) = dtiC * (-Ii_na - Ii_k - Ii_l - I_syn_ii - I_syn_ei - I_syn_ei_r - I_syn_ii_r + Ii - Ii_pulse + Ii_dyn(:,T));%
        Khi(:,rk) = dtiphi*(alpha_h .*(1-hih)- beta_h.*hih);
        Kni(:,rk) = dtiphi*(alpha_n .*(1-nih)- beta_n.*nih);
        Ksi(:,rk) = dt * (ialpha_syn * Fi .* (1-sih) - ibeta_syn * sih);
        
        KVe(:,rk) = dteC * (-Ie_na - Ie_nap - Ie_kdr - Ie_ka - Ie_l - I_syn_ee - I_syn_ie - I_syn_ee_r - I_syn_ie_r + Ie - Ie_pulse + Ie_dyn(:,T)); % - Ie_kslow(:)
        Kbe(:,rk) = dt * (be_inf - beh) / 15;
        Khe(:,rk) = dt * (he_inf - heh) ./ the;
        Kne(:,rk) = dt * (ne_inf - neh) ./ tne;
        Kse(:,rk) = dt * (ealpha_syn * Fe .* (1-seh) - ebeta_syn * seh);
    end
    
    si(:,Tsnew) = si(:, Tsprev) + (Ksi(:,1)+2*Ksi(:,2)+2*Ksi(:,3)+Ksi(:,4))/6;
    se(:,Tsnew) = se(:, Tsprev) + (Kse(:,1)+2*Kse(:,2)+2*Kse(:,3)+Kse(:,4))/6;
    
    ni = ni + (Kni(:,1)+2*Kni(:,2)+2*Kni(:,3)+Kni(:,4))/6;
    hi = hi + (Khi(:,1)+2*Khi(:,2)+2*Khi(:,3)+Khi(:,4))/6;
    
    be = be + (Kbe(:,1)+2*Kbe(:,2)+2*Kbe(:,3)+Kbe(:,4))/6;
    ne = ne + (Kne(:,1)+2*Kne(:,2)+2*Kne(:,3)+Kne(:,4))/6;
    he = he + (Khe(:,1)+2*Khe(:,2)+2*Khe(:,3)+Khe(:,4))/6;
    
    %noise
    noisei = -siIlambda + 2*siIlambda*rand(Ni*reg,1);
    noisee = -seIlambda + 2*seIlambda*rand(Ne*reg,1);
    
    Vi(:,Tvnew) = Vi(:,Tvprev)+(KVi(:,1)+2*KVi(:,2)+2*KVi(:,3)+KVi(:,4))/6 +sdt*noisei;
    Ve(:,Tvnew) = Ve(:,Tvprev)+(KVe(:,1)+2*KVe(:,2)+2*KVe(:,3)+KVe(:,4))/6 +sdt*noisee;
    
    
    %spikefinding!
    tempspikese{T} = find(Ve(:,Tvprev) <= thresholde & Ve(:,Tvnew) >= thresholde);
    tempspikesi{T} = find(Vi(:,Tvprev) <= thresholdi & Vi(:,Tvnew) >= thresholdi);
    
    spiketimese{T} = T*ones(size(tempspikese{T}));
    spiketimesi{T} = T*ones(size(tempspikesi{T}));

end

idata = cell2mat(tempspikesi);
edata = cell2mat(tempspikese);
spikesi = [mod(idata, Ni), (cell2mat(spiketimesi)-1)*dt, floor((idata-1)./Ni)+1];
spikese = [mod(edata, Ne), (cell2mat(spiketimese)-1)*dt, floor((edata-1)./Ne)+1];

end

