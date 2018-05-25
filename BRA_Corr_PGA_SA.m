function [rho] = BRA_Corr_PGA_SA(T)
% Correlation between PGA and horizontal spectral acceleration
%
%
% Implements the correlation equations developed by Bradley B. A. and
% published as "Empirical correlation of PGA, spectral accelerations and
% spectrum intensities from active shallowcrustal earthquakes" 
% (2011, Earthquake Engineering & Structural Dynamics pages 1707--1721)
%
%
if (0.2<=T) && (T<=10)
an=interp1([0.2;10],[1;0.97],T);
bn=interp1([0.2;10],[0.895;0.25],T);
cn=interp1([0.2;10],[0.06;0.80],T);
dn=interp1([0.2;10],[1.6;0.80],T);
end

if (0.2<=T) && (T<=10)
rho=real((an+bn)/2-(an-bn)/2*atanh(dn*log(T/cn)));
else
rho=interp1([0;0.2],[1;0.9173],T); 
end
end

