function [rho, p] = nancorr(sig1, sig2)

sel = ~isnan(sig1) & ~isnan(sig2);
[rho, p] = corr(sig1(sel), sig2(sel));

end