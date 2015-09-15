function getDishData_script
%% load g_clone and m_clone objects
load g_clone;
load m_clone;
g_lin=lineage(g_clone);
m_lin=lineage(m_clone);%% Analyses diffusion distances between cells in a dish and outputs dataframe for analysis in R
%% Analysis by R functions in METS package: cor.cif, Grandom.cif, com.risk (timereg package)
header={'Event_Time' 'Start_Time' 'Cause' 'Growth_Factor' 'Clone' 'Progenitor' 'Generation'};
g_listmode=ShowData(g_lin,'Listmode');
m_listmode=ShowData(m_lin,'Listmode');
  % rows are ...
    % 1. clone number
    % 2. cell ID (track)
    % 3. cell generation
    % 4. birth phenotype         
    % 5. death phenotype
    % 6. number of phenotype transitions per lifetime
    % 7. transition age
    % 8. lifetime 
    % 9. stop reason None=0(right-censored), Division, Apoptosis, Lost
    %10. BirthAge

outlier = 1e4; %microns
max_time = 6;  %days

%% plot positions of GFP positive cells
[g_dish] = PlotTrajectories( g_clone, g_listmode, outlier, max_time );
[m_dish] = PlotTrajectories( m_clone, m_listmode, outlier, max_time );

%% find all branches with GFP transition
g_GFP_ndx=~isnan(g_listmode.wavelength(7,:));
m_GFP_ndx=~isnan(m_listmode.wavelength(7,:));
% look up GFP transition intervals
[g_GFP_time,g_GFP_StartTime]=getTransitionInterval(g_listmode.wavelength,g_GFP_ndx,'GFP_on');
[m_GFP_time,m_GFP_StartTime]=getTransitionInterval(m_listmode.wavelength,m_GFP_ndx,'GFP_on');
%% look up right censored intervals
g_RC_ndx=(g_listmode.wavelength(9,:)==0)&~...
    (getDescendants(g_listmode.wavelength,g_GFP_ndx,0)); % right censored
g_RC_ndx=g_RC_ndx|...
    ((g_listmode.wavelength(9,:)==3)&~...
    (getDescendants(g_listmode.wavelength,g_GFP_ndx,3))); % lost as well
m_RC_ndx=(m_listmode.wavelength(9,:)==0)&~...
    (getDescendants(m_listmode.wavelength,m_GFP_ndx,0)); % right censored
m_RC_ndx=m_RC_ndx|...
    ((m_listmode.wavelength(9,:)==3)&~...
    (getDescendants(m_listmode.wavelength,m_GFP_ndx,3))); % lost as well
[g_RC_time,g_RC_StartTime]=getTransitionInterval(g_listmode.wavelength,g_RC_ndx,'Right_Censored');
[m_RC_time,m_RC_StartTime]=getTransitionInterval(m_listmode.wavelength,m_RC_ndx,'Right_Censored');
%% look up apoptotic intervals
g_apop_ndx=(g_listmode.wavelength(9,:)==2)&~...
    (getDescendants(g_listmode.wavelength,g_GFP_ndx,2)); % apoptosis
m_apop_ndx=(m_listmode.wavelength(9,:)==2)&~...
    (getDescendants(m_listmode.wavelength,m_GFP_ndx,2)); % right censored
[g_apop_time,g_apop_StartTime]=getTransitionInterval(g_listmode.wavelength,g_apop_ndx,'Death');
[m_apop_time,m_apop_StartTime]=getTransitionInterval(m_listmode.wavelength,m_apop_ndx,'Death');
%% assemble fields
Event_Time=[g_GFP_time(:);m_GFP_time(:);g_RC_time(:);m_RC_time(:);g_apop_time(:);...
    m_apop_time(:)];
Start_Time=[g_GFP_StartTime(:);m_GFP_StartTime(:);g_RC_StartTime(:);m_RC_StartTime(:);...
    g_apop_StartTime(:);...
    m_apop_StartTime(:)];
Cause=[ones(size(g_GFP_time(:)));ones(size(m_GFP_time(:)));zeros(size(g_RC_time(:)));...
    zeros(size(m_RC_time(:)));2*ones(size(g_apop_time(:)));2*ones(size(m_apop_time(:)))];
Clone=[g_listmode.wavelength(1,g_GFP_ndx) m_listmode.wavelength(1,m_GFP_ndx) ...
    g_listmode.wavelength(1,g_RC_ndx) m_listmode.wavelength(1,m_RC_ndx)...
    g_listmode.wavelength(1,g_apop_ndx) m_listmode.wavelength(1,m_apop_ndx)];
Clone=Clone';
Progenitor=[g_listmode.wavelength(2,g_GFP_ndx) m_listmode.wavelength(2,m_GFP_ndx) ...
    g_listmode.wavelength(2,g_RC_ndx) m_listmode.wavelength(2,m_RC_ndx)...
    g_listmode.wavelength(2,g_apop_ndx) m_listmode.wavelength(2,m_apop_ndx)];
Progenitor=Progenitor';
Generation=[g_listmode.wavelength(3,g_GFP_ndx) m_listmode.wavelength(3,m_GFP_ndx) ...
    g_listmode.wavelength(3,g_RC_ndx) m_listmode.wavelength(3,m_RC_ndx)...
    g_listmode.wavelength(3,g_apop_ndx) m_listmode.wavelength(3,m_apop_ndx)];
Generation=Generation';
Growth_Factor=[ones(size(g_GFP_time(:)));2*ones(size(m_GFP_time(:)));ones(size(g_RC_time(:)));...
    2*ones(size(m_RC_time(:)));ones(size(g_apop_time(:)));2*ones(size(m_apop_time(:)))];
%% Save as text file for R called "GFPregression"
%'Event_Time' 'Cause' 'Growth_Factor' 'Clone' 'Progenitor' 'Generation'
% get rid of NaNs
b=~isnan(Event_Time);
Event_Time=Event_Time(b);
Cause=Cause(b);
Growth_Factor=Growth_Factor(b);
Clone=Clone(b);
Progenitor=Progenitor(b);
Generation=Generation(b);
Start_Time=Start_Time(b);
SavePrompt='Save competing hazards regression data';
Matlab2CSV4R(header,SavePrompt,Event_Time,Start_Time,Cause,Growth_Factor,Clone,Progenitor,Generation);