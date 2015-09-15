function getDishData_script
%% load g_clone and m_clone objects
load g_clone;
load m_clone;
g_lin=lineage(g_clone);
m_lin=lineage(m_clone);%% Analyses diffusion distances between cells in a dish and outputs dataframe for analysis in R