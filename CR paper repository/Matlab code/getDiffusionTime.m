function [diff_t,del_x,LR]=getDiffusionTime(Pairs,T,CloneID,D,clones,n_dim,r)
%% Estimated expected diffusion time for a molecule with diffusion coefficient D
% Pairs are the progeny IDs that are compared, m x 2 array
% T is the event times, m x 2 array
% Calculates the likelihood ratio (LR) given a n_dim dimensional diffusion model
% with diffusion coeffcient D.
% The likelihood ratio is calculated by assuming that at r (one cell radius in microns) 
% the likelihood is unity at t=r/sqrt(D)
% diff_t is the cell separation (del_x) divided by the sqrt of the diffusion
% coefficient. 


n=size(Pairs,1);
del_x=zeros(n,1);del_t=zeros(n,1);LR=zeros(n,1);diff_t=zeros(n,1);
% is T relative to birth time or absolute time
if T(1)>datenum(1900,1,1)
    AbsoluteTime=true;
else
    AbsoluteTime=false;
end
for i=1:n
    Progeny=Pairs(i,:);
    trackA=Progeny2Track(clones{CloneID},Progeny(1));
    trackB=Progeny2Track(clones{CloneID},Progeny(2));
    trackA=clones{CloneID}.track{trackA};
    trackB=clones{CloneID}.track{trackB};
    %look up position when event occurs
    if AbsoluteTime
        ndx=find(trackA.X(1,:)>=T(i,1),1,'first');
        PosA=[trackA.X(2,ndx),trackA.Y(2,ndx)];
        ndx=find(trackB.X(1,:)>T(i,2),1,'first');
        PosB=[trackB.X(2,ndx),trackB.Y(2,ndx)];
    else
        ndx=find(trackA.X(1,:)>(T(i,1)+trackA.BirthTime),1,'first');
        PosA=[trackA.X(2,ndx),trackA.Y(2,ndx)];
        ndx=find(trackB.X(1,:)>(T(i,2)+trackB.BirthTime),1,'first');
        PosB=[trackB.X(2,ndx),trackB.Y(2,ndx)];
    end
    % for apoptosis or right censoring
    if isempty(PosA)
        pos=[trackA.X(2,:);trackA.Y(2,:)];
        b=(pos-repmat(pos(:,1),1,size(pos,2)))>500; % check for spurious x,y coordinates.
        b=~sum(b,1);
        pos=pos(:,b);
        if isempty(pos)
            PosA=[NaN,NaN];
        else
            PosA=[pos(1,end) pos(2,end)];
        end
    end
    if isempty(PosB)
        pos=[trackB.X(2,:);trackB.Y(2,:)];
        b=(pos-repmat(pos(:,1),1,size(pos,2)))>500; % check for spurious x,y coordinates.
        b=~sum(b,1);
        pos=pos(:,b);
        if isempty(pos)
            PosB=[NaN,NaN];
        else
            PosB=[pos(1,end) pos(2,end)];
        end
    end
    del_x(i)=norm(PosA-PosB); % in microns
    del_t(i)=abs(T(i,1)-T(i,2)); 
    diff_t(i)=(del_x(i)*1e-6)^2/D; % in seconds
    % convert to SI units
    if del_t(i)==0
        del_t(i)=1/24/60; %set to 1 minute
    end
    LR(i)=Density((PosA(1)-PosB(1))*1e-6,(PosA(2)-PosB(2))*1e-6,del_t(i)*24*60*60,D,n_dim);
    %normalise so likelihood ratio is greater than one inside cell at t=r/sqrt(D) 
    L_at_r=Density(r*1e-6/sqrt(2),r*1e-6/sqrt(2),(r*1e-6)^2/D,D,n_dim);
    LR(i)=LR(i)/L_at_r;
    if isnan(LR(i))
        disp('LR is NaN');
    end
 end

function Track=Progeny2Track(clone,Progeny)
n=length(clone.track);
found=false;
i=1;
while ~found&&(i<=n)
    if clone.track{i}.TrackNum==Progeny
        Track=i;
        found=true;
    else
        i=i+1;
    end
end

function G=Density(x,y,t,D,n)
% response to a dirac delta impulse for diffusion with infinite boundary
% conditions i.e. the green function is
% G(x,x',t)=exp(-|x-x'|^2/(4*D*t)) / (4*pi*D*t)^(n/2)
% for 2 dimensions
% assume x' = (0,0)
G=exp(-(x^2+y^2)/4/D/t)/(4*pi*D*t).^(n/2);





