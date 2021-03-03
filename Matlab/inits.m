clear all; clc;

% Initialization
COBRASOLVER = 'ibm_cplex'; %glpk | mosek | ibm_cplex | gurobi 'pdco'
%edit(fullfile(userpath,'startup.m'))
%setenv('ILOG_CPLEX_PATH','D:\Program Files\IBM\ILOG\CPLEX_Studio_Community1210')
% Add paths
addpath(genpath('datasets'));
addpath('methods');
addpath(genpath('utils'));

changeCobraSolver(COBRASOLVER, 'all');

clear COBRASOLVER ans;

