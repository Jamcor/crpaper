function cumlambda=getCSCumHazard(T,C)
%% get the Cause specific cumulative hazard where T is the time of
% failure and C is the cause, 0 reserved for censored data
X=[T(:) C(:)];
X=sortrows(X,1);
X=[0 0;X];
NumOfCauses=max(C);
ni=zeros(size(X,1),1);
d=zeros(length(ni),NumOfCauses);
for i=1:length(ni)
    ni(i)=risk(X(i,1),X(:,1));
    for j=1:NumOfCauses
        d(i,j)=fail(X(i,1),X(:,1),j,X(:,2));
    end
end
lambda=d./repmat(ni,1,NumOfCauses);
cumlambda=cumsum(lambda,1);
cumlambda=[X(:,1) cumlambda];
cumlambda=cumlambda(1:end-1,:);

    



function r=risk(t,T)
%% calculates risk set at time t using the vector time series T
r=sum(T>t);

function f=fail(t,T,c,C)
%% calculates the specific number of failures at time t using the 
% vector T (time series) and C (cause)
f=sum((T==t)&(C==c));
