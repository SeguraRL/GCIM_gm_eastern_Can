function [ rho ] = BRA_Corr_PGV_SA(T)
% Correlation between PGV and horizontal spectral acceleration
%
%
% Implements the correlation equations developed by Bradley B. A. and
% published as "Empirical correlations between peak ground velocity and 
% spectrum-based intensity measures" 
% (2012, Earthquake Spectra pages 17--35)
%
%
if (0.1<=T) && (T<0.75)
an=interp1([0.1;0.75],[0.73;0.54],T);
bn=interp1([0.1;0.75],[0.54;0.81],T);
cn=interp1([0.1;0.75],[0.045;0.28],T);
dn=interp1([0.1;0.75],[1.8;1.5],T);
elseif (0.75<=T) && (T<2.5)
an=interp1([0.75;2.5],[0.54;0.8],T);
bn=interp1([0.75;2.5],[0.81;0.76],T);
cn=interp1([0.75;2.5],[0.28;1.1],T);
dn=interp1([0.75;2.5],[1.5;3.0],T);
elseif (2.5<=T) && (T<=10)
an=interp1([2.5;10],[0.8;0.76],T);
bn=interp1([2.5;10],[0.76;0.70],T);
cn=interp1([2.5;10],[1.1;5.0],T);
dn=interp1([2.5;10],[3.0;3.2],T);
end

if (0.1<=T) && (T<=10)
rho=real((an+bn)/2-((an-bn)/2)*atanh(dn*log(T/cn)));
else
rho=interp1([0;0.1],[0.733;0.5534],T);
end
end


