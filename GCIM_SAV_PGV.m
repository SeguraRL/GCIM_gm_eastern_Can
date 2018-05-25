function [meanReq, covReq] = GCIM(T_star,T_h,T_v, M,Rcd,epsilon,VS30,PGA_1100);
%% Target event
% This function estimates the Mean and the Covariance matrix of vertical
% and horizontal SA at different periods and PGV conditioned on Sa_h(T1). The
% generalized conditional intensity measure approach proposed by Bradley
% (2012) is used. 

[median_sa_T_star, sigma_sa_T_star] = ATAD_2013_SA(M, T_star, Rcd);
lnSa1=median_sa_T_star;

%% Ground motion prediction Eq. for SA - Unconditional

for i=1:length(T_h)
   [median_sa(i), sigma_sa(i)] = ATAD_2013(M, T_h(i), Rcd); % median of Sa and standard deviation of lnSa
    rho_sa(i) = baker_jayaram_correlation(T_star, T_h(i)); % correlation coefficient with T*
end

%% Ground motion prediction Eq. for Sa vertical - Unconditional

for k=1:length(T_v)
   [m_VH sigma_VH] = GUL_2011_VH(M, T_v(k), Rcd, VS30,PGA_1100);
   [m_sav(k), s_sav(k)] = ATAD_2013_SA(M, T_v(k), Rcd); % median of Sa and standard deviation of lnSa
    rho_sav(k) = BKR_SAH_SAV(T_star,  T_v(k)); % correlation coefficient with T*
end
median_sav=m_sav.*m_VH;
sigma_sav=s_sav.*sigma_VH;

%% Ground motion prediction Eq. for PGV - Unconditional

T_pgv=999;
[median_pgv sigma_pgv]=ATAD_2013_pgv(M, T_pgv, Rcd);%PGV for the target scenario (cm/s)
rho_pgv=BRA_Corr_PGV_SA(T_star);% correlation coefficient with T*

%% PARAMETERIZATION OF THE TARGET RESPONSE SPECTRUM DISTRIBUTION 
%Mean conditioned on SA(T*)
mu_lnsa_c=log(median_sa)+rho_sa.*epsilon.*sigma_sa;
mu_lnsav_c=log(median_sav)+rho_sav.*epsilon.*sigma_sav;
mu_lnpgv_c=log(median_pgv)+rho_pgv.*epsilon.*sigma_pgv;

Mu=[mu_lnsa_c mu_lnsav_c mu_lnpgv_c]';

meanReq=Mu;
h=length(T_h);
v=length(T_v);

% Covariance matrix unconditioned 

for i=1:length(T_h)
for j=1:length(T_h)
for k=1:length(T_v)
for z=1:length(T_v)
rho_sa_unc(i,j)=baker_jayaram_correlation(T_h(i), T_h(j));
rho_sav_unc(k,z)=BKR_SAV_SAV(T_v(z),T_v(k));
rho_sav_sa_unc(i,k)=BKR_SAH_SAV(T_h(i),T_v(k));
rho_pgv_sa_unc(i)=BRA_Corr_PGV_SA(T_h(i));
rho_pgv_sav_unc(k)=BRA_Corr_PGV_SA(1.5).*BKR_SAH_SAV(1.5,T_v(k));

end
end
end
end

Sg=[sigma_sa sigma_sav sigma_pgv]';
D_sg=diag(Sg.^2);
rho_b=zeros(h+v+1);
rho_b(1:h,1:h)=rho_sa_unc;
rho_b(1:h,(h+1):(h+v))=rho_sav_sa_unc;
rho_b((h+1):(h+v),1:h)=rho_sav_sa_unc';
rho_b((h+1):(h+v),(h+1):(h+v))=rho_sav_unc;
rho_b((h+v+1),1:h)=rho_pgv_sa_unc;
rho_b(1:h,(h+v+1))=rho_pgv_sa_unc';
rho_b((h+v+1),(h+1):(h+v))=rho_pgv_sav_unc;
rho_b((h+1):(h+v),(h+v+1))=rho_pgv_sav_unc';
rho_b((h+v+1),(h+v+1))=1;

for u=1:(h+v+1)
for w=1:(h+v+1)    
    S_0(u,w)=rho_b(u,w)*Sg(u,1)*Sg(w,1);
end
end

% Covariance matrix between IMs and SA(T*) 
rho_cond_sa=[rho_sa';rho_sav';rho_pgv];
S_1=rho_cond_sa.*Sg.*sigma_sa_T_star;
 
% Covariance matrix conditioned on SA(T*)
SS=S_0-(1/(sigma_sa_T_star^2))*S_1*S_1';

% The conditional standard deviation
sigma_cond=Sg.*sqrt(1-rho_cond_sa.^2);

covReq_0=SS;
 
% The conditional correlation matrix

for u=1:(h+v+1)
for w=1:(h+v+1)   
    rho_cond(u,w)=(rho_b(u,w)-rho_cond_sa(u,1)*rho_cond_sa(w,1))/(sqrt(1-rho_cond_sa(u,1)^2)*sqrt(1-rho_cond_sa(w,1)^2));
end
end

e=eig(covReq_0);

% Determine the cholesky decomp of the corr matrix
    [L_b,p]=chol(covReq_0,'lower');
    
% check that corr matrix is PD and then get L
   
    if p~=0 %i.e. matrix is not PD
        options.keepDiag=true;
        %options.doDykstra = true;
        %options.maxit=500;
        [output]=NearPD(covReq_0,options);
        if output.converged~=true
            fprintf('Error: NearPD.m failed to return PD correlation matrix \n');
            return;
        end
        
% determine the Frobenius norm of the nearPD matrix
        covReq_PD=output.X;
        [L]=chol(covReq_PD,'lower');
        normF=output.normF;
    end
    
        
    if any(e<0)==1
        covReq=covReq_PD;
    else
        covReq=covReq_0;
    end

