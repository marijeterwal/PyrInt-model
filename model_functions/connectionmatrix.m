function [conee_out, conei_out, conii_out, conie_out, coneer_out, coneir_out, coniir_out, conier_out] = ...
    connectionmatrix(sett, outputind)

Ni = sett.Ni;
Ne = sett.Ne;
reg = sett.Nregions;

%% Local:
s = RandStream('mcg16807','Seed',sett.seed+15);
RandStream.setGlobalStream(s);

if sett.connectiontype == 1; %all to all
    pl = [1, 1, 1, 1];
elseif sett.connectiontype == 2; %  no connections
    pl = [0,0,0,0];
elseif sett.connectiontype == 3; % random #2
    pl = [sett.Msynii, sett.Msynie, sett.Msynee, sett.Msynei];
end

conii = zeros(Ni,Ni,reg); conii(rand(Ni,Ni,reg) <= pl(1)) = 1;
conie = zeros(Ne,Ni,reg); conie(rand(Ne,Ni,reg) <= pl(2)) = 1;
conee = zeros(Ne,Ne,reg); conee(rand(Ne,Ne,reg) <= pl(3)) = 1;
conei = zeros(Ni,Ne,reg); conei(rand(Ni,Ne,reg) <= pl(4)) = 1;

%% Global:
s = RandStream('mcg16807','Seed',sett.seed+16);
RandStream.setGlobalStream(s);

if reg >1;
    
    if length(sett.Msynee_rtor(1,:)) < reg || length(sett.Msynee_rtor(:,1)) < reg || ...
            length(sett.Msynei_rtor(1,:)) < reg || length(sett.Msynei_rtor(:,1)) < reg;
        error('Msyn definition are invalid: should be a reg*reg matrix.');
    end
    
    if length(sett.g_syn_ee_r(1,:)) < reg || length(sett.g_syn_ee_r(:,1)) < reg || ...
            length(sett.g_syn_ei_r(1,:)) < reg || length(sett.g_syn_ei_r(:,1)) < reg;
        error('Gsyn definition are invalid: should be a reg*reg matrix.');
    end
    
    ncomb = 2*nchoosek(reg,2);
    comb = [combnk(1:reg,2); fliplr(combnk(1:reg,2))];
    
    conee_r = zeros(Ne,Ne,ncomb);
    conei_r = zeros(Ni,Ne,ncomb);
    conii_r = zeros(Ni,Ni,ncomb);
    conie_r = zeros(Ne,Ni,ncomb);
    
    for nc = 1:ncomb
        reg1 = comb(nc,1); % from
        reg2 = comb(nc,2); % to
        
        if sett.connectiontype == 1; %all to all
            pg = [1,1,1,1];
        elseif sett.connectiontype == 2; %  no connections
            pg = [0,0,0,0];
        elseif sett.connectiontype == 3; % random #2
            
            % test diagonals
            if reg1 == reg2 && sum(pg) > 0; error('Connection probabilities are invalid. Diagonals should be zero.'); end
            if reg1 == reg2 && sett.g_syn_ee_r(reg2,reg1) == 0; error('Connection strengths are invalid. Diagonals should be zero.'); end
            if reg1 == reg2 && sett.g_syn_ei_r(reg2,reg1) == 0; error('Connection strengths are invalid. Diagonals should be zero.'); end
            if reg1 == reg2 && sett.g_syn_ii_r(reg2,reg1) == 0; error('Connection strengths are invalid. Diagonals should be zero.'); end
            if reg1 == reg2 && sett.g_syn_ie_r(reg2,reg1) == 0; error('Connection strengths are invalid. Diagonals should be zero.'); end
            
            pg = [sett.Msynee_rtor(reg2, reg1), sett.Msynei_rtor(reg2, reg1)...
                sett.Msynii_rtor(reg2, reg1), sett.Msynie_rtor(reg2, reg1)];
        end
        
        conee_rtemp = zeros(Ne,Ne,1); conee_rtemp(rand(Ne,Ne,1) <= pg(1)) = 1;
        conei_rtemp = zeros(Ni,Ne,1); conei_rtemp(rand(Ni,Ne,1) <= pg(2)) = 1;
        conii_rtemp = zeros(Ni,Ni,1); conii_rtemp(rand(Ni,Ni,1) <= pg(3)) = 1;
        conie_rtemp = zeros(Ne,Ni,1); conie_rtemp(rand(Ne,Ni,1) <= pg(4)) = 1;
        conee_r(:,:,nc) = conee_rtemp;
        conei_r(:,:,nc) = conei_rtemp;
        conii_r(:,:,nc) = conii_rtemp;
        conie_r(:,:,nc) = conie_rtemp;
    end
end

%% Number of connections

Nii_l = mean(sum(conii, 2));
Nie_l = mean(sum(conie, 2));
Nee_l = mean(sum(conee, 2));
Nei_l = mean(sum(conei, 2));
for r1 = 1:reg
    if Nii_l(:,:,r1) == 0 || isnan(Nii_l(:,:,r1)); Nii_l(:,:,r1) = 1; end
    if Nie_l(:,:,r1) == 0 || isnan(Nie_l(:,:,r1)); Nie_l(:,:,r1) = 1; end
    if Nee_l(:,:,r1) == 0 || isnan(Nee_l(:,:,r1)); Nee_l(:,:,r1) = 1; end
    if Nei_l(:,:,r1) == 0 || isnan(Nei_l(:,:,r1)); Nei_l(:,:,r1) = 1; end
end

if reg > 1
    Nee_r = mean(sum(conee_r, 2));
    Nei_r = mean(sum(conei_r, 2));
    Nii_r = mean(sum(conii_r, 2));
    Nie_r = mean(sum(conie_r, 2));
    for nc = 1:ncomb
        if Nee_r(:,:,nc) == 0 || isnan(Nee_r(:,:,nc)); Nee_r(:,:,nc) = 1; end
        if Nei_r(:,:,nc) == 0 || isnan(Nei_r(:,:,nc)); Nei_r(:,:,nc) = 1; end
        if Nii_r(:,:,nc) == 0 || isnan(Nii_r(:,:,nc)); Nii_r(:,:,nc) = 1; end
        if Nie_r(:,:,nc) == 0 || isnan(Nie_r(:,:,nc)); Nie_r(:,:,nc) = 1; end
    end
end

%% combine conee/conei and conee_r/conei_r

if outputind == 2; %plotting
    conii_out = zeros(Ni*reg, Ni*reg);
    conie_out = zeros(Ne*reg, Ni*reg);
    conee_out = zeros(Ne*reg, Ne*reg);
    conei_out = zeros(Ni*reg, Ne*reg);
    
    for r = 1:reg
        conii_out((r-1)*Ni+1:r*Ni,(r-1)*Ni+1:r*Ni) = conii(:,:,r);
        conie_out((r-1)*Ne+1:r*Ne,(r-1)*Ni+1:r*Ni) = conie(:,:,r);
        conee_out((r-1)*Ne+1:r*Ne,(r-1)*Ne+1:r*Ne) = conee(:,:,r);
        conei_out((r-1)*Ni+1:r*Ni,(r-1)*Ne+1:r*Ne) = conei(:,:,r);
    end
    
    if reg >1;
        for nc = 1:ncomb;
            reg1 = comb(nc,1); % from
            reg2 = comb(nc,2); % to
            conee_out((reg2-1)*Ne+1:reg2*Ne, (reg1-1)*Ne+1:reg1*Ne) = conee_r(:,:,nc);
            conei_out((reg2-1)*Ni+1:reg2*Ni, (reg1-1)*Ne+1:reg1*Ne) = conei_r(:,:,nc);
            %
        end
    end
    coneer_out = 0;
    coneir_out = 0;
    coniir_out = 0;
    conier_out = 0;
    
elseif outputind == 1; %calc
    %% g_syn / N *con (to speed up calculations in calc function)
    
    g_syn_ee_l = sett.g_syn_ee;
    g_syn_ei_l = sett.g_syn_ei;
    g_syn_ii_l = sett.g_syn_ii;
    g_syn_ie_l = sett.g_syn_ie;
    g_syn_ee_r = sett.g_syn_ee_r;
    g_syn_ei_r = sett.g_syn_ei_r;
    g_syn_ii_r = sett.g_syn_ii_r;
    g_syn_ie_r = sett.g_syn_ie_r;
    
    conii_out = zeros(Ni*reg, Ni*reg);
    conie_out = zeros(Ne*reg, Ni*reg);
    conee_out = zeros(Ne*reg, Ne*reg);
    conei_out = zeros(Ni*reg, Ne*reg);
    coneer_out = zeros(Ne*reg, Ne*reg); % matrix with to*from
    coneir_out = zeros(Ni*reg, Ne*reg); % matrix with to*from
    coniir_out = zeros(Ni*reg, Ni*reg); % matrix with to*from
    conier_out = zeros(Ne*reg, Ni*reg); % matrix with to*from
    
    for r = 1:reg %to
        conii_out((r-1)*Ni+1:r*Ni,(r-1)*Ni+1:r*Ni) = (g_syn_ii_l ./ Nii_l(:,:,r)) .* conii(:,:,r);
        conie_out((r-1)*Ne+1:r*Ne,(r-1)*Ni+1:r*Ni) = (g_syn_ie_l ./ Nie_l(:,:,r)) .* conie(:,:,r);
        conee_out((r-1)*Ne+1:r*Ne,(r-1)*Ne+1:r*Ne) = (g_syn_ee_l ./ Nee_l(:,:,r)) .* conee(:,:,r);
        conei_out((r-1)*Ni+1:r*Ni,(r-1)*Ne+1:r*Ne) = (g_syn_ei_l ./ Nei_l(:,:,r)) .* conei(:,:,r);
    end
    if reg >1;
        for nc = 1:ncomb
            reg1 = comb(nc,1); % from
            reg2 = comb(nc,2); % to
            
            coneer_out((reg2-1)*Ne+1:reg2*Ne, (reg1-1)*Ne+1:reg1*Ne) = (g_syn_ee_r(reg2,reg1)  ./ Nee_r(:,:,nc)) .* conee_r(:,:,nc);
            coneir_out((reg2-1)*Ni+1:reg2*Ni, (reg1-1)*Ne+1:reg1*Ne) = (g_syn_ei_r(reg2,reg1)  ./ Nei_r(:,:,nc)) .* conei_r(:,:,nc);
            coniir_out((reg2-1)*Ni+1:reg2*Ni, (reg1-1)*Ni+1:reg1*Ni) = (g_syn_ii_r(reg2,reg1)  ./ Nii_r(:,:,nc)) .* conii_r(:,:,nc);
            conier_out((reg2-1)*Ne+1:reg2*Ne, (reg1-1)*Ni+1:reg1*Ni) = (g_syn_ie_r(reg2,reg1)  ./ Nie_r(:,:,nc)) .* conie_r(:,:,nc);
            
        end
    end
    
end

conii_out = sparse(conii_out);
conie_out = sparse(conie_out);
conee_out = sparse(conee_out);
conei_out = sparse(conei_out);
coneer_out = sparse(coneer_out);
coneir_out = sparse(coneir_out);
coniir_out = sparse(coniir_out);
conier_out = sparse(conier_out);
end