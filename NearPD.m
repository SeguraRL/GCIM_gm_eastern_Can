function [output]=NearPD(x,options)
% by Bradley B.A, 2012
% Purpose: This function applies the method of Higham (2002) to find the
%'nearest' PD matrix.
%
%reference: Higham, N. J., 2002. Computing the nearest correlation matrix
%—a problem from finance, IMA Journal of Numerical Analysis,  22, 329-343.
%
%Input variables:
%Rho - correlation matrix which is not SPD
%options - a structure with the following inputs
%       .maxits
%
%Output variables:
%Rho_mod - SPD correlation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start of core code
debug=false;

%defaults
options=setDefaults(options);
keepDiag=options.keepDiag;
doDykstra = options.doDykstra;
eigTol = options.eigTol;
convTol = options.convTol;
maxit =  options.maxit;

n=size(x,1);
D_S = zeros(n,n);

if keepDiag==true
    for i=1:n
        diagX0(i) = x(i,i);
    end
end

%make copy of original matrix
X = x;

%Set iteration, convergence criteria
iter = 0;
converged = false;
conv = Inf('double');
%



% [L,p]=chol(Rho);
% if p==0
%     SPD = 'yes';
% else
%     SPD = 'no';
% end

%print to screen only if debug mode
% if debug==true;  fprintf('SPD: %s \n',SPD); end;

% it=0;
while (iter<maxit)&(converged==false)     %(p>0)&(it<maxits)
    Y = X;
%     it=it+1;

    %Dykstra correction
    if (doDykstra)
        R = Y - D_S;
    end
    
   %project onto PSD matrices  X_k  =  P_S (R_k)
    if (doDykstra)
        [V,D] = eig(R);
    else 
        [V,D] = eig(Y);
    end
    eigen = diag(D);
    
    d = real(eigen);
    Q = V;
	%Get the maximum eigenvalue
	eigMax=-conv;
    for i=1:n
        if d(i)>eigMax
            eigMax = d(i);
        end
    end
	%compute the D_plus diagonal matricies
    for i=1:n
        d_plus = max(d(i),eigTol*eigMax);
        D_plus(i,i)=d_plus;
    end
	        
    X = (Q*D_plus)*(Q');
	        
    %Update Dykstra correction
    if (doDykstra)
        D_S = X-R;
    end

    %project onto symmetric and possibly 'given diag' matrices:
    if (keepDiag)
        for i=1:n
            X(i,i)=diagX0(i);
        end
    end

    %update convergence and iteration values
    conv = norm(Y-X,inf) / norm(Y,inf);
    iter = iter + 1;

    %check convergence criteria
    if (conv <= convTol)
        converged=true;
    end
end
	
%Set solution local variables as globals
output.X = X;
output.conv = conv;
output.normF = norm(x-X,'fro');
output.iter = iter;
output.eigVals = d;
output.converged=true;
    
%     %compute the eigenvalues
%     d=eig(Rho);
%     [V,D]=eig(Rho);
%     
%     if debug==true; fprintf('Setting the negative eig vals to +eps \n'); end;
%     
%     for d=1:size(D,1); 
%         if D(d,d)<0; D(d,d)=eps; end;
%     end
%     Rho_mod=V*D*inv(V);
%     
%     if debug==true; 
%         fprintf('Rho modified \n');
%     
%         for d=1:size(D,1);
%             for dd=1:size(D,1)
%                 fprintf(' %5.3f',Rho_mod(d,dd));
%             end
%             fprintf('\n');
%         end
% 
%         fprintf('Rho original \n');
%         for d=1:size(D,1);
%             for dd=1:size(D,1)
%                 fprintf(' %5.3f',Rho(d,dd));
%             end
%             fprintf('\n');
%         end
%     end
%     
%     %now check that Rho_mod is PD
%     [L,p]=chol(Rho_mod);
%     if p==0; 
%         SPD = 'yes';
%     else
%         SPD = 'no';
%         Rho = Rho_mod;
%     end
%     %print to screen
%     if debug==true; fprintf('SPD Rho_mod: %s \n',SPD); end

%end of function    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = setDefaults(options)
%Sets defaults

if isfield(options,'keepDiag')~=1
    options.keepDiag=false;
end

if isfield(options,'doDykstra')~=1
    options.doDykstra = true;
end

if isfield(options,'eigTol')~=1
    options.eigTol=1.e-6;
end

if isfield(options,'convTol')~=1
    options.convTol=1.e-7;
end

if isfield(options,'maxit')~=1
    options.maxit=100;
end

%end of function    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%