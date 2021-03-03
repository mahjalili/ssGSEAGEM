clear all; clc;

% Initialization
COBRASOLVER = 'ibm_cplex'; % glpk | mosek | ibm_cplex | gurobi | pdco
changeCobraSolver(COBRASOLVER, 'all');
% Add paths
addpath(genpath('datasets'));
addpath('methods');
addpath(genpath('utils'));

clear COBRASOLVER ans;
