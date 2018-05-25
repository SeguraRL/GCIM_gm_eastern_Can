function [ Sa Sigma ] = ATAD_2013_pgv(M, T, Rcd)
% by Rocio Segura, 18/10/2016
%
% Implements GMPE developed by Atkinson, G.M. and Adams J. and published as
% “Ground motion prediction equations for application to the 2015 Canadian
% national seismic hazard maps” 
% (2013, Canadian Journal of Civil Engineering, pages 988--998). 

% Originally we should consider the GMMs for eastern Canada from Atkinson
% and Boore (2006), Silva (2002) and Pezeshk et al. (2011) but the latter
% does not provide the values of PGV.

%No site parameters are needed. The GMPE was developed for hard-rock site
% with Vs30 >= 2000 m/s (NEHRP site class A) only.
[Sa_2, sigma_2]=AB_2006A_ENA(M, T, Rcd);
[Sa_3, sigma_3]=SLV_2002_DOUBLE_ENA(M, T, Rcd);
[Sa_4, sigma_4]=SLV_2002_SINGLE_ENA(M, T, Rcd);
A=[Sa_2  Sa_3  Sa_4]';
B=[sigma_2  sigma_3  sigma_4]';
Sa=mean(A);
Sigma=mean(B);

end

