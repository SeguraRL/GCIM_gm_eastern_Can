function [ Sa Sigma ] = ATAD_2013_SA(M, T, Rcd)
%
% Implements GMPE developed by Atkinson, G.M. and Adams J. and published as
% “Ground motion prediction equations for application to the 2015 Canadian
% national seismic hazard maps” 
% (2013, Canadian Journal of Civil Engineering, pages 988--998). 
%
[Sa_1, sigma_1]=PZT_2011_ENA(M, T, Rcd);
[Sa_2, sigma_2]=AB_2006A_ENA(M, T, Rcd);
[Sa_3, sigma_3]=SLV_2002_DOUBLE_ENA(M, T, Rcd);
[Sa_4, sigma_4]=SLV_2002_SINGLE_ENA(M, T, Rcd);
A=[Sa_1  Sa_2  Sa_3  Sa_4]';
B=[sigma_1  sigma_2  sigma_3  sigma_4]';
Sa=mean(A);
Sigma=mean(B);

end

