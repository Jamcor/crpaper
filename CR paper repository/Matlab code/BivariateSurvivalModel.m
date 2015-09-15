function [ X,Y,theta, rho, pval,rho_est,S,C ] = BivariateSurvivalModel( tau,mit_h,apop_h,sigma,N,cens_h)
%Bivariate survival model to test competing risks regression analysis
%   MonteCarlo simulation

%   The joint Survival model is a Gamma frailty model where theta is gamma
%   distributed i.e. cumumative hazard (ch) for t1 and t2 are
%    ch1(t1)=theta*mit_h*(t1>tau)(t1-tau)  
%    ch1(t2)=theta*apop_h*t2

%   mit is the mitosis hazard, tau is the S phase lag
%   apop is the apoptosis hazard
%   alpha=beta=1/sigma^2 are the parameters of the Gamma distibution
%   cens_h is the independent censoring hazard
%   N is the number of simulated pairs

%   X and Y are the pre-competition/censoring mitosis and apoptosis times
%   theta is the frailty
%   rho is Kendall's tau
%   rho_est is calculated using the Clayton model i.e.
%   rho_est=sigma^3/(sigma^2+2)
%   S is the survivor time
%   C is the cause (0 censored, 1 mitosis, 2 apoptosis)

theta=random('Gamma',1/(sigma^2),sigma^2,[N,1]); %mean 1, variance sigma^2

if abs(mean(theta)-1)>0.1
    theta=theta-mean(theta)+1;
    display('Sigma too low for random number generator');
end
% simulate mitosis
X=random('exp',1./(theta*mit_h))+tau;
Y=random('exp',1./(theta*apop_h));
[rho,pval]=corr(X,Y,'type','Kendall');
psi=sigma^2;
rho_est=psi/(psi+2);

% now compete apoptosis and mitosis
S=min(X,Y);
censoring=random('exp',1/cens_h,N,1);
C(S==X)=1;
C(S==Y)=2;
C(S>censoring)=0;
b=S>censoring;
S(b)=censoring(b);
if abs(rho_est-rho)/rho>0.2
    display(['Inaccurate simulation, rho=' num2str(rho) ', rho_est=' num2str(rho_est)]);
end
end

