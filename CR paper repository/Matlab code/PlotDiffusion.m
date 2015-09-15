function PlotDiffusion(D,r,n,d,t)
% r and d in microns, t in seconds, D in SI units
% D=1e-10; r=10; d=[2 5 10 20 40 100 200]; t=[0.01 20000]; n=3;
r=r*1e-6; %convert to metres
T=logspace(log10(t(1)),log10(t(2)),1000);
close all
d=d*1e-6; % in meteres
Y=r/sqrt(n)*ones(1,n);
L_at_r=Density(Y,r.^2/D,D,n);
Y=d./sqrt(n);
Y=repmat(Y(:),[1,n]);

for i=1:size(Y,1)
    LR(:,i)=Density(Y(i,:),T(:),D,n)/L_at_r;
    str{i}=[num2str(d(i)*1e6) ' um'];
end
plot(T,LR);
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'YLim',[1e-5,1e3]);
set(gca,'XLim',[1e-3,3e3]);
legend(str);
xlabel('seconds')
ylabel('Likelihood ratio');







function G=Density(X,t,D,n)
% response to a dirac delta impulse for diffusion with infinite boundary
% conditions i.e. the green function is
% G(X,t)=exp(-|X|^2/(4*D*t)) / (4*pi*D*t)^(n/2)
% X is a m -by- n array of position vectors
% for n dimensions
G=exp(-sum(X.^2,2)./4./D./t)./((4.*pi.*D.*t).^(n/2));