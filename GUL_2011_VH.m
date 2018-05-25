function [m_VH sigma_VH] = GUL_2011_VH(M, T, Rcd, VS30,PGA_1100)
% by Rocio Segura, 26/10/2016
% 
% Implements the GMMs for the vertical-to-horizontal spectral acceleration
% (V/H) ratio published as "Site-Specific Design Spectra for Vertical
% Ground Motion." (2011, Earthquake Spectra, pages 1023--1047)
%
%
%**************************************************************************
%--------------------------INPUT VARIABLES---------------------------------
% M      = moment magnitude
% T      = period (s); 
% Rcd    = closest distance to fault (km)
% VS30   = average shear-wave velocity over the top 30 m (m/s)
% PGA_1100 = is the median peak horizontal acceleration (g) for VS30 = 1100 m/s
%--------------------------OUTPUT VARIABLES--------------------------------
% VH    = median V/H ratio
% sigma = standard deviation of the V/H ratio
%**************************************************************************

%------Period-independent constants for the median V/H ratio --------------
% Table of coefficients were provided by Table 2.
c1=6.75;
c4=10;
a3=0.0147;
a4=0.0334;
a5=-0.034;
n=1.18;
co=1.88;

%----------------------COEFFICIENTS FOR ROCK SITES ------------------------
% Table of coefficients were provided by Table 6.
period = [0.01	0.02	0.029	0.04	0.05	0.075	0.1	0.15  0.2	0.26	0.3	 0.4	0.5  0.75	1 1.5  2	3	4	5	7.5	 10]';
  
cr = [865.1	-1.186	0.14	-0.16	-0.105	0	0.003	-1.23	0.45	0.33	0.23	0.15
      865.1	-1.219	0.14	-0.16	-0.105	0	0.003	-1.268	0.45	0.33	0.23	0.15
      898.6	-1.269	0.335	-0.185	-0.14	0	0.003	-1.366	0.45	0.33	0.23	0.15
      994.5	-1.308	0.562	-0.238	-0.16	0	0.003	-1.457	0.45	0.341	0.23	0.15
      1053.5 -1.346	0.72	-0.275	-0.136	0	-0.001	-1.533	0.45	0.351	0.23	0.15
      1085.7 -1.471	0.552	-0.24	-0.019	0	-0.007	-1.706	0.45	0.37	0.23	0.15
      1032.5 -1.624	0.214	-0.169	0	0.017	-0.01	-1.831	0.45	0.384	0.23	0.15
      877.6	-1.931	-0.262	-0.069	0	0.04	-0.008	-2.114	0.45	0.403	0.23	0.15
      748.2	-2.188	-0.6	0.002	0	0.057	-0.003	-2.362	0.45	0.416	0.23	0.15
      639	-2.412	-0.769	0.023	0	0.072	0.001	-2.527	0.45	0.429	0.23	0.15
      587.1	-2.518	-0.861	0.034	0	0.08	0.006	-2.598	0.45	0.436	0.23	0.15
      503	-2.657	-1.045	0.057	0	0.097	0.015	-2.685	0.45	0.449	0.23	0.15
      456.6	-2.669	-1.189	0.075	0	0.11	0.022	-2.657	0.45	0.46	0.23	0.15
      410.5	-2.401	-1.25	0.09	0	0.133	0.022	-2.265	0.45	0.479	0.237	0.15
      400	-1.955	-1.209	0.09	0	0.15	0.022	-1.685	0.45	0.492	0.266	0.15
      400	-1.025	-1.152	0.09	0.029	0.15	0.022	-0.57	0.45	0.511	0.307	0.15
      400	-0.299	-1.111	0.09	0.05	0.15	0.022	0.25	0.532	0.52	0.337	0.15
      400	0	-1.054	0.09	0.079	0.15	0.022	0.46	0.648	0.52	0.378	0.213
      400	0	-1.014	0.09	0.1	0.15	0.022	0.46	0.7	0.52	0.407	0.258
      400	0	-1	0.09	0.1	0.15	0.022	0.46	0.7	0.52	0.43	0.292
      400	0	-1	0.09	0.1	0.15	0.022	0.46	0.7	0.52	0.471	0.355
      400	0	-1	0.09	0.1	0.15	0.022	0.46	0.7	0.52	0.5	0.4];   
  

%**************************************************************************
% FIND C COEFFICIENTS BY INTERPOLATING BETWEEN PERDIODS
%**************************************************************************
ilow = min(find(period<=T));
T_low = period(ilow);
ihigh = max(find(period>=T));
T_high = period(ihigh);
% if given period equals a period in the table, then no need to interpolate
if ihigh==ilow
        c = cr(ihigh,:);
% otherwise, interpolate between coeffients
else
        c_high = cr(ihigh,:);
        c_low = cr(ilow,:);
    for i=1:length(c_high)
        c(i) = interp1([T_low T_high], [c_low(i) c_high(i)], T);
    end
end

Vlin=c(1);
b=c(2);
a1=c(3);
a2=c(4);
a6=c(5);
a7=c(6);
a8=c(7);
a10=c(8);
s1=c(9);
s2=c(10);
s3=c(11);
s4=c(12);

%**************************************************************************
% COMPUTE f1
%**************************************************************************

R=(Rcd^2+c4^2)^0.5;
if M<=c1
    f1 = a1+a4*(M-c1)+a8*(8.5-M)^2+(a2+a3*(M-c1))*log(R);
else
    f1 = a1+a5*(M-c1)+a8*(8.5-M)^2+(a2+a3*(M-c1))*log(R);
end

%**************************************************************************
% COMPUTE f5
%**************************************************************************

if T<=0.5
    V1=1500;
elseif 0.5<T<=1
    V1=exp(8.0-0.795*log(T/0.21));
elseif 1<T<2
    V1=exp(6.76-0.297*log(T));
else
    V1=700;
end

if VS30<V1
    VS30_star=VS30;
else
    VS30_star=VS30;
end
    
if VS30_star<Vlin
    f5=a10*log(VS30_star/Vlin)-b*log(PGA_1100+co)+b*log(PGA_1100+co*(VS30_star/Vlin)^n);
else
    f5=a10*log(VS30_star/Vlin)-(b*n)*log(VS30_star/Vlin);
end
    
%**************************************************************************
% COMPUTE F_NM OR F_RV
% FNM is the flag for normal faulting earthquakes (1 for normal earthquakes defined by
% rake angles between ?60 and ?120 degrees, 0 otherwise) and FRV is the flag for reverse-
% faulting earthquakes (1 for reverse or reverse/oblique earthquakes defined by rake angles between
% 30 and 150 degrees, 0 otherwise).
%**************************************************************************
F_NM=1;
F_RV=0;

%**************************************************************************
% Logarithm of the V/H ratio
%**************************************************************************
ln_VH = f1+a6*F_RV+a7*F_NM+f5;

%**************************************************************************
% CONVERT SPECTRAL ACCELERATION FROM LOG SCALE 
%**************************************************************************
m_VH = (exp(ln_VH));

%**************************************************************************
% DEFINE SIGMA TOTAL
%**************************************************************************
if M<5
    sigma_o=s1;
elseif 5<=M<=7
    sigma_o=s1+(s2-s1)*0.5*(M-5);
else
    sigma_o=s2;
end

if M<5
    tao_o=s3;
elseif 5<=M<=7
    tao_o=s3+(s4-s3)*0.5*(M-5);
else
    tao_o=s4;
end

sigma_VH=(sigma_o^2+tao_o^2)^0.5;
