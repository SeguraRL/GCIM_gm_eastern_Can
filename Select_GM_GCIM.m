% This code is used to select conditional (structure- and site-
% specific) ground motions. The target means and covariances are
% obtained corresponding to a pre-defined target scenario earthquake, and
% are obtained based on the GCIM method
%
% Base code:
% Nirmal Jayaram, Ting Lin, Jack W. Baker
% Department of Civil and Environmental Engineering
% Stanford University
% Last Updated: 21 July 2010
%
% Modified by: Rocio Segura
% Université de Shebrooke
% Department of Civil and Environmental Engineering
% Last Updated: 05 March 2018
%
%% Reference manuscripts:
%
% B.A. Bradley (2010). A generalized conditional intensity measure approach
% and holistic ground-motion selection, Earthquake Engineering and
% Structural Dynamics.
%
% B.A. Bradley (2012). A ground motion selection algorithm based on the
% generalized conditional intensity measure approach, Soil Dynamics and 
% Earthquake Engineering.
%
% Baker, J. W. and Lee, C. (2016). “An Improved Algorithm for Selecting
% Ground Motions to Match a Conditional Spectrum.” 
% Journal of Earthquake Engineering 
%
% N. Jayaram, T. Lin and and Baker, J. W. (2011). A computationally
% efficient ground-motion selection algorithm for matching a target
% response spectrum mean and variance, Earthquake Spectra.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT VARIABLES
% finalRecords      : Record numbers of selected records
% finalScaleFactors : Scale factors
%
% The final cell in this m file shows how to plot the selected spectra
% using this information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATABASE
% The NGA-WEST2 database is used. For an alternate database, the minimum 
% information to be provided includes the pseudo-acceleration spectra of 
% the available, ground motions, periods at which the spectra are defined,
% and other information % required to compute means and variances using 
% ground-motion models.
% This cell can be modified by the user if desired.
%
% Variable definitions
% saKnown   : (N*P matrix)
%             This is a matrix of Sa values at different periods (P) for
%             available ground-motion time histories (N).
% perKnown  : The set of P periods.
% nGM       : Number of ground motions to be selected
% T1        : Period at which spectra should be scaled and matched.
%             Usually, the structure's fundamental period.
% isScaled  : Should spectra be scaled before matching (1 -YES, 0-NO).
% maxScale  : Maximum allowable scale factor.
% weights   : [Weight for error in mean, Weight for error in standard
%             deviation] e.g., [1.0,1.0] - equal weight for both errors.
% nLoop     : This is a meta-variable of the algorithm. The greedy
%             improvement technique does several passes over the available
%             data set and tries to improve the selected records. This
%             variable denotes the number of passes. Recommended value: 2.
% penalty   : If a penalty needs to be applied to avoid selecting spectra
%             that have spectral acceleration values beyond 3 sigma at any
%             of the periods, set a value here. Use 0 otherwise.
% notAllowed: List of record numbers that should not be considered for
%             selection. Use [] to consider all available records. This can
%             be used to prevent certain magnitude-distance-Vs30 records
%             from being selected, if desired. (see example below)
% seedValue : For repeatability. For a particular seedValue not equal to
%             zero, the code will output the same set of ground motions.
%             The set will change when the seedValue changes. If set to
%             zero, the code randomizes the algorithm and different sets of
%             ground motions (satisfying the target mean and variance) are
%             generated each time.
% outputFile: File name of the output file
%
% If a database other than the NGA database is used, also define the
% following variables:
%
% magnitude        : Magnitude of all the records
% distance_closest : Closest distance for all the records
% soil_Vs30        : Soil Vs30 values corresponding to all the records
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs begin here
clear; close all; clc
load NGA_W2_meta_data % Information from the NGA database
load NGA_W2_PGV % Vel in cm/s
SaKnown    =Sa_RotD50(logical(PGV>0),:);
PGV_new    = PGV(logical(PGV>0));
Sa_vert   = abs(Sa_vert);
Sa_vert_m   = Sa_vert(logical(PGV>0),:);
perKnown   = Periods;
nGM        = 30;
T1         = 0.25;
Sa_T1      = 0.2; %(g) Target intensity level
isScaled   = 1;
maxScale   = 5;
weights    = [1.0 1.0];
nLoop      = 0;
penalty    = 0;
notAllowed = 1;
seedValue  = 1;
outputFile = 'Output_File.dat';
% NOTE: MORE user input required in the next cell

% Limiting the records to be considered using the 'notAllowed' variable

% Only hard rock
recInvalid = find(soil_Vs30<550 | soil_Vs30>2000);
notAllowed = [notAllowed; recInvalid];

% Limits on magnitude
recInvalid = find(magnitude<5.4 | magnitude>7.5);
notAllowed = [notAllowed; recInvalid];

% Limits on distance
recInvalid = find(closest_D<0 | closest_D>200);%
notAllowed = [notAllowed; recInvalid];

% MORE user input in the next cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determination of target mean and covariances
% Variable definitions

% M_bar     : Magnitude of the target scenario earthquake
% R_bar     : Distance corresponding to the target scenario earthquake
% eps_bar   : Epsilon value for the target scenario
% Vs30      : Average shear wave velocity in the top 30m of the soil, used
%             to model local site effects (m/s)
% Ztor      : Depth to the top of coseismic rupture (km)
% delta     : Average dip of the rupture place (degree)
% lambda    : Rake angle (degree)
% Zvs       : Depth to the 2.5 km/s shear-wave velocity horizon (km)
% arb       : 1 for arbitrary component sigma
%             0 for average component sigma
% PerTgt    : Periods at which the target spectrum needs to be computed
% showPlots : 1 to display plots, 0 otherwise
% useVar    : 0 to set target variance to zero, 1 to compute target
%             variance using ground-motion model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs begin here

M_bar     = 6.01;
R_bar     = 54.18;
eps_bar   = 0.78; 
Vs30      = 1500;
PGA_1100  = 0.302;
Ztor      = 0;
delta     = 90;
lambda    = 180;
Zvs       = 2;
arb       = 0;
T_v       = logspace(log10(0.2*T1),log10(2*T1),20); % periods of interest in the vertical direction (0.2T1-2T1)
T_h       = logspace(log10(0.2*T1),log10(2*T1),20);

%Modify T_h to include T1
if ~any(T_h == T1)
T_h=[T_h(T_h<T1) T1 T_h(T_h>T1)];
end 
h=length(T_h);
v=length(T_v);
PerTgt    = T_h;
PerTgt_total    = [T_h T_v];
w_h       = 0.60/h+zeros(1,h);
w_v       = 0.20/v+zeros(1,v);
w_pgv       = 0.20;
w_IM      = [w_h w_v w_pgv];
showPlots = 1;
useVar    = 1;

Rrup   = R_bar; % Can be modified by the user
Rjb    = R_bar; % Can be modified by the user

% User inputs end here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GCIM Mean and Covariance
[meanReq, covReq] = GCIM_SAV_PGV(T1,T_h,T_v, M_bar,R_bar,eps_bar,Vs30,PGA_1100);
[idx1 idx1] = min(abs(T_h-T1)); %index of closest value
 lnSa1 = meanReq(idx1,1);
 meanReq=meanReq';
 Ss=sqrt(diag(covReq));
%% Simulate response spectra using Monte Carlo Simulation

% Sets of response spectra are simulated and the best set (in terms of
% matching means, variances and skewness is chosen as the seed). 

nTrials = nGM; 
% Setting initial seed for simulation
if seedValue ~= 0
    randn('seed',seedValue);
else
    randn('seed',sum(100*clock));
end
devTotalSim = zeros(nTrials,1);
for j=1:nTrials
    gmCell{j} = zeros(nGM,length(meanReq));
    for i=1:nGM
    gmCell{j}(i,:) = exp(mvnrnd(meanReq,covReq));
    end
    devMeanSim = mean(log(gmCell{j})) - meanReq;
    devSkewSim = skewness(log(gmCell{j}),1);
    devSigSim = std(log(gmCell{j})) - sqrt(diag(covReq))';
    devTotalSim(j) = weights(1) * sum(devMeanSim.^2) + weights(2) * sum(devSigSim.^2)+ 0.1 * (weights(1)+weights(2)) * sum(devSkewSim.^2);
end
[tmp recUse] = min(abs(devTotalSim));
gm = gmCell{recUse};

if showPlots == 1
   
    
%% Plot simulated response spectra
    gf_sim=figure;
    loglog(PerTgt, exp(meanReq(1,1:h)), '-r', 'linewidth', 3)
    hold on
    loglog(PerTgt, exp(meanReq(1,1:h) + 1.96*sqrt(diag(covReq(1:h,1:h)))'), '--r', 'linewidth', 3)
    for i=1:nGM
        loglog(PerTgt,(gm(:,1:h))','k');
    end
    loglog(PerTgt, exp(meanReq(1,1:h) - 1.96*sqrt(diag(covReq(1:h,1:h)))'), '--r', 'linewidth', 3)
    axis([min(PerTgt) max(PerTgt) 1e-2 5])
    xlabel('T (s)')
    ylabel('S_a (g)')
    legend('Median response spectrum','2.5 and 97.5 percentile response spectra','Response spectra of simulated ground motions')
    set(gca, 'fontsize', 14, 'FontName','Euclid');
    
%% Plot target and simulated means
    mean_sim=figure;
    loglog(PerTgt,exp(meanReq(1,1:h)))
    hold on
    loglog(PerTgt,exp(mean(log(gm(:,1:h)))),'--')
    axis([min(PerTgt) max(PerTgt) 1e-2 5])
    xlabel('T (s)')
    ylabel('Median S_a (g)')
    legend('exp(Target mean lnS_a)','exp(Mean of simulated lnS_a)')
    title('Target and sample exponential logarithmic means (i.e., medians)')
    set(gca, 'fontsize', 14, 'FontName','Euclid');
    
%% Plot target and simulated standard deviations
    std_sim=figure;
    semilogx(PerTgt,sqrt(diag(covReq(1:h,1:h)))')
    hold on
    semilogx(PerTgt,std(log(gm(:,1:h))),'--')
    axis([min(PerTgt) max(PerTgt) 0 1])
    xlabel('T (s)')
    ylabel('Standard deviation of lnS_a')
    legend('Target standard deviation of lnS_a','Standard deviation of simulated lnS_a')
    title('Target and sample logarithmic standard deviations')
    set(gca, 'fontsize', 14, 'FontName','Euclid');

    
end

%% Arrange the available spectra in a usable format and check for invalid input

% Match periods (known periods and periods for error computations)
recPer_h = zeros(length(T_h),1);
for i=1:length(T_h)
    [tmp recPer_h(i)] = min(abs(perKnown - T_h(i)));
end

recPer_v = zeros(length(T_v),1);
for i=1:length(T_v)
    [tmp_v recPer_v(i)] = min(abs(perKnown - T_v(i)));
end

% Check for invalid input
sampleBig_h = SaKnown(:,recPer_h);
sampleBig_v = Sa_vert_m(:,recPer_v);
sampleBig_pgv =PGV_new;
sampleBig = [sampleBig_h sampleBig_v sampleBig_pgv];
if (any(any(isnan(sampleBig))))
    error ('NaNs found in input response spectra')
end

% Processing available spectra
sampleBig = log(sampleBig);
nBig = size(sampleBig,1);

%% Find best matches to the simulated spectra from ground-motion database

recID = zeros(nGM,1);
sampleSmall = [];
finalScaleFac = ones(nGM,1);
for i = 1:nGM
    err = zeros(nBig,1);
    
    scaleFac = ones(nBig,1);
    for j=1:nBig
    if (isScaled == 1)
            
    if exp(sampleBig(j,idx1)) == 0
       scaleFac(j) = -1;
       err(j) = 10e6;
    else
       scaleFac(j) = exp(lnSa1)/exp(sampleBig(j,idx1));
    if (scaleFac(j) > maxScale | soil_Vs30(j)==-1 | any(notAllowed==j))
        err(j) = 10e6;
    else
        err(j) = w_IM*(((log(exp(sampleBig(j,:))*scaleFac(j)) - log(gm(i,:)))./Ss').^2)';
    end
    end
    else
    if (soil_Vs30(j)==-1 | any(notAllowed==j))
       err(j) = 10e6;
    else
       err(j) = w_IM*(((log(exp(sampleBig(j,:))*scaleFac(j)) - log(gm(i,:)))./Ss').^2)';
    if err(j) == inf
       err(j) = 10e6;
    end
    end
    end
    if (any(recID == j))
        err(j) = 10e6;
    end
    end
    [tmp recID(i)] = min(err);
    if tmp >= 10e6
        display('Warning: Possible problem with simulated spectrum. No good matches found');
        display(recID(i));
    end
    if (isScaled == 1)
        finalScaleFac(i) = scaleFac(recID(i));
    else
        finalScaleFac(i) = 1;
    end
    sampleSmall = [sampleSmall;log(exp(sampleBig(recID(i),:))*scaleFac(recID(i)))];
    
end


%% Greedy subset modification procedure

display('Please wait...This algorithm takes a few minutes depending on the number of records to be selected');

for k=1:nLoop % Number of passes
    
for i=1:nGM % Selects nSelect ground motions
        
    display([num2str(round(((k-1)*nGM + i-1)/(nLoop*nGM)*100)) '% done']);
        
    minDev = 10e6;
        
    sampleSmall(i,:) = [];
    recID(i,:) = [];
        
% Try to add a new spectra to the subset list
for j=1:nBig
            
if isScaled == 1
if exp(sampleBig(j,idx1)) == 0
   scaleFac(j) = 10e6;
else
   scaleFac(j) = exp(lnSa1)/exp(sampleBig(j,idx1));
end
   sampleSmall = [sampleSmall;sampleBig(j,:)+log(scaleFac(j))];
else
   sampleSmall = [sampleSmall;sampleBig(j,:)];
   scaleFac(j) = 1;
end
% Compute deviations from target
devMean = mean(sampleSmall) - meanReq;
devSkew = skewness(sampleSmall,1);
devSig = std(sampleSmall) - sqrt(diag(covReq))';
devTotal = weights(1) * sum(devMean.^2) + weights(2) * sum(devSig.^2);
            
% Penalize bad spectra (set penalty to zero if this is not required)
for m=1:size(sampleSmall,1)
    devTotal = devTotal + sum(abs(exp(sampleSmall(m,:))>exp(meanReq+3*sqrt(diag(covReq))'))) * penalty;
end
if (scaleFac(j) > maxScale | soil_Vs30(j)==-1 | any(notAllowed==j))
   devTotal = devTotal + 10e6;
end
            
% Should cause improvement and record should not be repeated
if (devTotal < minDev && ~any(recID == j))
    minID = j;
    minDev = devTotal;
end
sampleSmall = sampleSmall(1:end-1,:);
end
        
% Add new element in the right slot
if isScaled == 1
   finalScaleFac(i) = scaleFac(minID);
else
   finalScaleFac(i) = 1;
end
   sampleSmall = [sampleSmall(1:i-1,:);sampleBig(minID,:)+log(scaleFac(minID));sampleSmall(i:end,:)];
   recID = [recID(1:i-1);minID;recID(i:end)];
end
end

display('100% done');

% Output information
finalRecords = recID;
finalScaleFactors = finalScaleFac;

%% Spectral Plots

if (showPlots)
        
    
    % Plot at all periods
    gf_selec=figure;
    loglog(PerTgt, exp(meanReq(1,1:h)), 'b', 'linewidth', 3)
    hold on
    loglog(PerTgt, exp(meanReq(1,1:h) + 1.96*sqrt(diag(covReq(1:h,1:h)))'), '--b', 'linewidth', 3)
    loglog(perKnown,SaKnown(finalRecords,:).*repmat(finalScaleFactors,1,size(SaKnown,2)),'color',[0.5 0.5 0.5]);
    loglog(PerTgt, exp(meanReq(1,1:h) - 1.96*sqrt(diag(covReq(1:h,1:h)))'), '--b', 'linewidth', 3)
    axis([min(PerTgt) max(PerTgt) 1e-2 5])
    xlabel('T (s)');
    ylabel('S_a (g)');
    legend('Median response spectrum','2.5 and 97.5 percentile response spectra','Response spectra of selected ground motions');
    set(gca, 'fontsize', 14, 'FontName','Euclid');

    
    % Plot spectra only at periods where error is minimized
    gf_selec=figure;
    loglog(PerTgt, exp(meanReq(1,1:h)), 'b', 'linewidth', 3)
    hold on
    loglog(PerTgt, exp(meanReq(1,1:h) + 1.96*sqrt(diag(covReq(1:h,1:h)))'), '--b', 'linewidth', 2.5)
    loglog(PerTgt,exp(sampleBig(finalRecords,1:h)).*repmat(finalScaleFactors,1,length(PerTgt)),'color',[0.5 0.5 0.5],'linewidth',1)
    loglog(PerTgt, exp(meanReq(1,1:h) - 1.96*sqrt(diag(covReq(1:h,1:h)))'), '--b', 'linewidth', 2.5)
    axis([min(PerTgt) max(PerTgt) 1e-2 5])
    xlabel('T (s)');
    ylabel('Sa (g)');
    axis([0.075 0.5 0.01 1])
    legend({'Mean response spectrum','2.5 and 97.5 percentile response spectra','Response spectra of selected GMTH'},'Box','off');
    set(gca, 'fontsize', 14, 'FontName','Euclid');
   
    
    % Sample and target means
    gf_selec_means=figure;
    loglog(PerTgt,exp(meanReq(1,1:h)),'k','linewidth',2)
    hold on
    loglog(PerTgt,exp(mean(sampleSmall(:,1:h))),'b--','linewidth',2)
    axis([min(PerTgt) max(PerTgt) 1e-2 5])
    xlabel('T (s)')
    ylabel('Median S_a (g)')
    legend('exp(Target mean lnS_a)','exp(Mean of selected lnS_a)')
    set(gca, 'fontsize', 14, 'FontName','Euclid');
    
    % Sample and target standard deviations
    gf_selec_dev=figure;
    semilogx(PerTgt,sqrt(diag(covReq(1:h,1:h)))','k','linewidth',2)
    hold on
    semilogx(PerTgt,std(sampleSmall(:,1:h)),'b--','linewidth',2)
    axis([min(PerTgt) max(PerTgt) 0 1])
    xlabel('T (s)')
    ylabel('Standard deviation of lnS_a')
    legend('Target standard deviation of lnS_a','Standard deviation of selected lnS_a','Location','northeast')
    set(gca, 'fontsize', 14, 'FontName','Euclid');
    

%% Selected unscaled records and horizontal scaled factor
sel_rec=finalRecords;
scale_fc=finalScaleFac;

%Selected vertical ground motions 
    
Sa_vert_selected=Sa_vert_m(sel_rec,:);
  
%Scale vertical ground motions with horizontal scale factor
for k=1:length(scale_fc);
    Sa_vert_scaled(k,:)=Sa_vert_selected(k,:).*scale_fc(k);
end
    
%Mean of selected vertical ground motions
for j=1:nGM
for w=1:length(perKnown)
if Sa_vert_scaled(j,w)<0
   Sa_vert_scaled(j,w)=NaN;
end
end
end
Mean_Sa_vert= nanmean(Sa_vert_scaled);

%% Kolmogorov–Smirnov (K–S) goodness-of-fit test for Sa_V(T1)

%Find Mu and sigma for an specific period
tmp = abs(T_v-T1);
[idx idx] = min(tmp); %index of closest value in the evrtical range
closest = T_v(idx); %closest value
mu_sav=meanReq(h+idx);
s_p_sav=sqrt(covReq(h+idx,h+idx));
T_s       = 0:0.001:(Sa_T1+0.1);
y_sav = cdf('Lognormal',T_s,mu_sav,s_p_sav);
targetCDF=makedist('Lognormal','mu',mu_sav,'sigma',s_p_sav);
[H,P,KSSTAT,CV] = kstest(gm(:,h+idx),'CDF',targetCDF,'Alpha',0.1); %just used to get CV nothing else



% Figure marginal IM and record selection
gf_CDF_sav=figure;
g=cdfplot(Sa_vert_scaled(:,logical(perKnown==T1)));
hold on
set(g,'color','b','linewidth', 1)
plot(T_s,y_sav,'r','linewidth', 2.5);
plot(T_s,y_sav+CV,'--r','linewidth', 1.75);
plot(T_s,y_sav-CV,'--r','linewidth', 1.75);
legend({'GMTH - GCIM','F_{Sa_V(T1)|Sa_H(T1)}','K-S bounds'},'Box','off','Location','southeast')
nx=xlabel('Sa_v (g)', 'FontSize', 14);
ny=ylabel('Cumulative Probability', 'FontSize', 14);
axis([0 0.5*(Sa_T1) 0 1])
set(gca, 'fontsize', 14, 'FontName','Euclid');
title('');

% Plot selected response spectra for the vertical component
gf_vert_mean=figure;
loglog(T_v, exp(meanReq(1,h+1:h+v)), '-r', 'linewidth', 2.5);
hold on
loglog(T_v, exp(meanReq(1,h+1:h+v) + 1.96*Ss(h+1:h+v,1)'), '--r', 'linewidth', 2.5)
plot(perKnown, Sa_vert_scaled(2:end,:), 'color',[0,0,0]+0.5, 'linewidth', .75);
loglog(T_v, exp(meanReq(1,h+1:h+v) - 1.96*Ss(h+1:h+v,1)'), '--r', 'linewidth', 2.5)
xlabel('T (s)')
ylabel('Sa_v (g)')
axis([0.3*T1 2*T1 0.001 1])
legend({'Median Sa_v(T)|Sa_h(T_1)','2.5 and 97.5 percentile Sa_v(T)|Sa_h(T_1)',...
        'Vertical Sa of selected GMTH'},'Box','off','Location','southeast')
set(gca, 'fontsize', 14, 'FontName','Euclid');

%% Kolmogorov–Smirnov (K–S) goodness-of-fit test for PGV
    
pgv_selected=PGV_new(sel_rec,:);
  
%Scale PGV with horizontal scale factor

pgv_scaled=pgv_selected.*scale_fc;

mu_pgv=meanReq(h+v+1);
s_p_pgv=sqrt(covReq(h+v+1,h+v+1));
T_s_pgv       = 0:.01:Sa_T1*50;
y_pgv = cdf('Lognormal',T_s_pgv,mu_pgv,s_p_pgv);
targetCDF=makedist('Lognormal','mu',mu_pgv,'sigma',s_p_pgv);
[H_pgv,P_pgv,KSSTAT_pgv,CV_pgv] = kstest(gm(:,h+v+1),'CDF',targetCDF,'Alpha',0.1);


% Figure marginal IM and record selection
gf_CDF_pgv=figure;
g=cdfplot(pgv_scaled);
hold on
set(g,'color','b','linewidth', 1)
plot(T_s_pgv,y_pgv,'r','linewidth', 2.5);
plot(T_s_pgv,y_pgv+CV_pgv,'--r','linewidth', 1.75);
plot(T_s_pgv,y_pgv-CV_pgv,'--r','linewidth', 1.75);
legend({'GMTH-GCIM','F_{PGV|Sa_H(T1)}','K-S bounds'},'Box','off','Location','southeast')
nx_1=xlabel('PGV(cm/s)', 'FontSize', 14);
ny_1=ylabel('Cumulative Probability', 'FontSize', 14);
axis([0 Sa_T1*40 0 1])
set(gca, 'fontsize', 14, 'FontName','Euclid');
title('');

end

%% Output data to file (best viewed with textpad)

%Records number in the global system

Rec_global=find(PGV>0);
finalRecords_mod=Rec_global(finalRecords,1);

fin = fopen(outputFile,'w');
fprintf(fin,'%s \t %s \t %s \t %s \t %s \t %s \t %s \n','Record Number','NGA Record Sequence Number','Scale Factor','File Name Dir. 1','File Name Dir. 2','URL Dir 1','URL Dir 2');
for i = 1 : length(finalRecords_mod)
    rec = finalRecords_mod(i);
    url1 = ['http://peer.berkeley.edu/nga_files/ath/' Filename_1{rec}(1:end-3) 'AT2'];
    url2 = ['http://peer.berkeley.edu/nga_files/ath/' Filename_2{rec}(1:end-3) 'AT2'];
    fprintf(fin,'%d \t %d \t %6.2f \t %s \t %s \t %s \t %s \n',i,rec,finalScaleFactors(i),Filename_1{rec},Filename_2{rec},url1,url2);
end

fclose(fin);