function [rho] = BKR_SAV_SAV(T1,T2)
% by Rocio Segura, 26/10/2016
% 
% Correlation between the  values of vertical spectral ground motion 
% component at two different periods.
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
% rho    = Correlation between vertical spectral acceleration at different
% periods
%**************************************************************************

T_max=max(T1,T2);
T_min=min(T1,T2);
rho=1-0.77*log(T_max/T_min)+0.315*(log(T_max/T_min))^1.4;
