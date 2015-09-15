
%% Load undifferentiated HSC clone file
load('D:\Cell tracking data\Jonathan Draper\Undiff hESC 4x 15min\-03040_-03918_4_f2\colony_1.mat');
%% Create sibling data.frame for R
% aim is to perform cross odds ratio using distance to edge of colony as a
% covariate for the marginal distribution of cell generation time
[Undiff_HSC_Sisters]=getSisters(clone);
save('D:\Cell tracking data\Jonathan Draper\Undiff hESC 4x 15min\SisterData','Undiff_HSC_Sisters');
