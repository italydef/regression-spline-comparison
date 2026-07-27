
function fit = bayes_pspline(x,y,nknots,nsubknots,degree,nIter,burn_in) ;


% This function fits a univariate smoothed regression spline using a P-spline basis
% uses MCMC to sample the parameters
% 
% INPUTS 
%
% REQUIRED:
% 
% x         = independent variables (n by 1)
% y         = vector of dependent variable (n by 1)
% nknots    = number of knots 
% nsubknots = number of subknots for penalty weights 
% 
% OPTIONAL:
%
% degree    = degree of p-spline (default=1)
% nIter     = number of iterations for MCMC (default=5000)
% burn_in   = burn in for the the MCMC (default = 1000)


% OUTPUTS
%   Returns a structure with the following fields -
%                             
% yhat        = estimates using the regression spline
% yhatlower   = lower 95% Bayesian limit for yhat
% yhatupper   = upper 95% Bayesian limit for yhat
% Sigma       = smoothing parameter at each knot i.e. variance of the  \Omega_y
% Omega_y     = values of \Omega_y generated at each iteration of MCMC (matrix with each row containing
%                                                                                the vector)
% knots
% subknots

% Calls: srslocal,gen_basis_pspline
% 
% Modified: 7/29/03  

if nargin<=4
    degree=1;
    nIter=5000;
    burn_in=1000;
end;

thining = 5;

%--------------------------------------------------------------------------------------

n=length(x);    % length of data vector
My =nknots;     % no of knots;
Ms = nsubknots; %no of knots for D matrix;
met = 5;        % number of MH replications within the MH algorithm
sample=(nIter-burn_in)/thining; %number of MCMC samples to be collected
%--------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Generation of starting values for MH step for \sigma^2(X) using code from %%%%%
%%%%%  Ruppert and Carroll (1999). The MATLAB codes can be downloaded from       %%%%%  
%%%%%  http://www.orie.cornell.edu/~davidr/  and has to stored in the same       %%%%%  
%%%%%  directory from which this function is called.                             %%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd srs
ruppert_est = srslocal(x,y,degree,My,Ms) ;
cd ..

knots1=ruppert_est.knots;
knots2=ruppert_est.subknots;

alpha_y = interp1(knots2,exp(ruppert_est.alpha)',knots1,'linear')' ;

sigmas = 1./(alpha_y);


%--------------------------------------------------------------------------------------

%%%%%  Generation of the P-spline basis matrix %%%%%


Zy=gen_pspline_basis(x,knots1,degree); %  For \Omegay_Y
Zs=gen_pspline_basis(knots1,knots2,degree); %  For \Omegay_s


%--------------------------------------------------------------------------------------


% Matrices in which values for parmaters generated at each MCMC iteration will be stored
% MUST FOR CONVERGENCE DIAGNOSTICS!!!!!!

omegay_mat= zeros(sample,1+degree+My);  % for \Omegay_Y
sigy_mat = zeros(sample,1);      % for sigma^2_Y
xi_mat = zeros(sample,1);        % for \Xi^2 
omegas_mat = zeros(sample,Ms+degree+1); % for \Omegay_s
sigmas_mat = zeros(sample,My);   % for \sigma^2(X)

%--------------------------------------------------------------------------------------

% Starting values for the MCMC sampler and values of hyperparameters

sig_alpha_y=100*ones(1,degree+1);    % for variance of fixed effects of first stage p-spline
sig_alpha_s=100*ones(1,degree+1);    % for variance of fixed effects of first stage p-spline

ay=10;by=0.01;      % Inverse Gamma Hyperparametes for sigma^2_Y
as=1000;bs=2;     % Inverse Gamma Hyperparametes for \Xi^2   

sigy=0.01;          % for sigma^2_Y
xi=0.01;            % for \Xi^2 

omegay=0.01*ones(My+degree+1,1);       % for \Omegay_Y
omegas=0.01*ones(Ms+degree+1,1);       % for \Omegay_s           
sig_u=0.01;                            % for sigma^2_u (note this is fixed in our calculations).

%--------------------------------------------------------------------------------------

lambdas=diag([sig_alpha_s xi*ones(1,Ms)]);  % Var-Cov matrix for \Omegay_s
lambday=diag([sig_alpha_y sigmas]);         % Var-Cov matrix for \Omegay_Y
ZytZy = (Zy' * Zy);                                                                      
ZstZs = (Zs' * Zs);     
flag=zeros(1,My);   % To monitor movement of MH step: 0 - No move ; 1 - Move


rho=-log(sigmas);

cd srs
ruppert_est2 = srslocal(knots1,rho',degree,Ms); 
cd ..

m = zeros(Ms+2,1);
 
%--------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                          START OF MCMC SAMPLER                             %%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;indx=0;indx1=0; %counters to be used inside the MCMC loop

while indx < nIter,
    
   count=count+1;
    %---------------for \Omegay_Y -------------;

  
        lambday=diag([sig_alpha_y sigmas]);
        zigmay= inv((ZytZy/sigy)+inv(lambday));
        mean_omegay=(zigmay*Zy'*y)/sigy; 
        omegay = (mvnrnd(mean_omegay,zigmay,1))';   
   
   
    %---------------for \sigma^2_Y -------------;
 
        temp=((y-Zy*omegay)'*(y-Zy*omegay))/2;  %';
        newby=1/(temp+(1/by));
        sigy = 1/(gamrnd((ay+(0.5*n)),newby));
       
 
   % ---------------for \Xi^2 -----------------;
     
        temp1=sum(((omegas(degree+2:Ms+degree+1))-m(degree+2:Ms+degree+1)).^2);
        newbs=1/((temp1/2)+(1/bs));
        xi = 1/(gamrnd((as+(0.5*Ms)),newbs));
        
   
 
    % ---------- for \Omegay_s -----------------;

        
        lambdas=diag([sig_alpha_s xi*ones(1,Ms)]);
        zigmas= inv((ZstZs/sig_u)+inv(lambdas));
        postmean_omegas=zigmas*(((Zs'*(-log(sigmas))')/sig_u) +inv(lambdas)*m);
        omegas = (mvnrnd(postmean_omegas,zigmas,1))'; 
        
        
    %---------- Metropolis-Hastings step withing Gibbs sampler -------------;
    %----------             for for \sigma^2(X)                -------------;
    
        betay=omegay(degree+2:My+degree+1);
        mean_log_sigmas=Zs*omegas;
        for i=1:My,
            old_sigmas=sigmas(i);
            old_eita = log(old_sigmas);
            old_likl= -(0.5)*(-old_eita+(((betay(i))^2)/exp(old_eita))+(((-old_eita-(mean_log_sigmas(i)))^2)/sig_u));
            indicator=1; % 1 for rejection 
              
            for j=1:met,
                new_eita=normrnd(old_eita,.0001);
                
                if(indicator==0)
                    old_likl=new_likl;;        
                end;
                
                new_likl=-(0.5)*(-new_eita+(((betay(i))^2)/exp(new_eita))+(((-new_eita-(mean_log_sigmas(i)))^2)/sig_u));
                a=new_likl-old_likl;
                u=log(unifrnd(0,1));
                
                if(u<a)
                    old_eita=new_eita;
                   indicator=0;
                 end;
       
            end; 
            
            old_sigmas=exp(old_eita);
                if(sigmas(i)~=old_sigmas)
                    flag(i)=flag(i)+1;
                end;
            sigmas(i)=old_sigmas;
        end;
        
             
        %------- Collecting post burn in samples ------;
        
        
        
        if count>burn_in, %start collecting samples            
            if(rem(count,thining)==0),
             indx1=indx1+1;
             omegay_mat(indx1,:)=omegay';
             sigy_mat(indx1,:)=sigy;
             xi_mat(indx1,:)=xi;
             omegas_mat(indx1,:)=omegas';
             sigmas_mat(indx1,:)=sigmas;
            end; 
        end;
        
          if rem(count,100)==0
            fprintf('Iteration: %d Samples : %d/%d  \n',count,indx1*5,sample*5);
          end
        indx = indx+1;
        
    end;  % End of MCMC sampler 

 
         
    
%--------------------------------------------------------------------------

yhat=mean(Zy*(omegay_mat'),2); % fitted values 

yhattemp=Zy*(omegay_mat');
yhatsorted=sort(yhattemp,2);

alpha=5; % alpha-level for Bayesian credible intervals

yhatlower=(prctile(yhatsorted',alpha))'; 
yhatupper=(prctile(yhatsorted',100-alpha))';

sigmastemp=sigmas_mat;
%sigmastemp(1:burn_in,:)=[]; %throwing out the the burn_in values
sigmas = mean(sigmastemp,1); 

%--------------------------------------------------------------------------------------


    fit=struct('x',x,...
              'y',y,...
              'yhat',yhat,...
              'yhatlower',yhatlower,...
              'yhatupper',yhatupper,...
              'Sigma',sigmas,...
              'knots',knots1,...
              'subknots',knots2,...
              'Omega_y',omegay_mat,...
              'Omega_s',omegas_mat,... 
              'Sigma_s',sigmas_mat);
                            