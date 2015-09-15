%% Frailty model script to test CR regression
tau=0.4;mit_h=2;apop_h=1/100;sigma=1e-4;N=1000;cens_h=1/100;
%% test the effects of apoptosis and mitosis hazard on the cumulative incidence estimator of Gray and Fine
Event_Time=[];
Cause=[];
Mitosis_Hazard=[];
Apoptosis_Hazard=[];
apop_h=[0,1];
for i=1:2        
    [ X,Y,theta, rho, pval,rho_est,S,C ] =...
        BivariateSurvivalModel( tau,mit_h,apop_h(i),sigma,N,cens_h);
    Mitosis_Hazard=[Mitosis_Hazard; mit_h*ones(N,1)];
    Apoptosis_Hazard=[Apoptosis_Hazard; apop_h(i)*ones(N,1)];
    Cause=[Cause;C(:)];
    Event_Time=[Event_Time;S(:)];      

end

header={'Event_Time' 'Cause' 'Mitosis_Hazard' 'Apoptosis_Hazard'};
Matlab2CSV4R(header,Event_Time,Cause,Mitosis_Hazard,Apoptosis_Hazard);