function [rho] = BKR_SAH_SAV(T1,T2)
% by Rocio Segura, 26/10/2016

% Correlation between the  values of vertical and
% horizontal spectral ground motion component at two different periods.
% 
% Implements the correlation equations developed by Baker and Cornell and
% published as "Correlation of response spectral values for multicomponent
% ground motions" (2006, Bulletin of the Seismological Society of America
% pages 215--227)
%
%**************************************************************************
%--------------------------INPUT VARIABLES---------------------------------
% T1,T2  = Periods of interest
%--------------------------OUTPUT VARIABLES--------------------------------
% rho    = Correlation between vertical and horizontal spectral acceleration 
% at different periods
%**************************************************************************

T_max=max(T1,T2);
T_min=min(T1,T2);
%Indicator function
if T_min<0.189
    I=1;
else
    I=0;
end
rho=(0.64+0.021*log(sqrt(T_min*T_max)))*(1-cos(0.5*pi-(log(T_max/T_min))*(0.29+0.094*I*(log(T_min/0.189)))));
