% import GMP lineage switching stats.txt via excell
%% plot cumulative incidence
clear all
load('GFPregressionStats');
h1=plot(cuminc(:,1),cuminc(:,2:5),'-');
hold on;
h2=plot(predcuminc_GFP(:,1),predcuminc_GFP(:,2:3),':',...
    predcuminc_apoptosis(:,1),predcuminc_apoptosis(:,2:3),':');
hold off;
legend([h1;h2],{'GCSF->GFP' 'MCSF->GFP' 'GCSF->apoptosis' 'MCSF->apoptosis',...
    'GCSF->GFP (reg)' 'MCSF->GFP (reg)' 'GCSF->apoptosis (reg)' 'MCSF->apoptosis (reg)'})
xlabel('days');
ylabel('Cumulative incidence');
set(gcf,'Name','Effect of GF on GFP expression')

%% calculate cumulative hazard
GCSF_CH=getCSCumHazard(GCSF_TS(:,1),GCSF_TS(:,2));
MCSF_CH=getCSCumHazard(MCSF_TS(:,1),MCSF_TS(:,2));
figure('Name','Cumulative Hazard');
plot(GCSF_CH(:,1),GCSF_CH(:,2:3),MCSF_CH(:,1),MCSF_CH(:,2:3));
legend('GCSF->GFP','GCSF->apoptosis','MCSF->GFP','MCSF->apoptosis');
xlabel('days');
ylabel('Cause specific cumulative hazard');