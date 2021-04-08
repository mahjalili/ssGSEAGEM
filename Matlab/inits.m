% Initialization
% A conceptual framework for transcriptional data integration into a genome-scale metabolic model.
%
% Author: Mahdi Jalili, 2021

clear all; clc;

COBRASOLVER = 'ibm_cplex'; % glpk | mosek | ibm_cplex | gurobi | pdco
changeCobraSolver(COBRASOLVER, 'all');
% Add paths
addpath(genpath('datasets'));
addpath('methods');
addpath(genpath('utils'));

clear COBRASOLVER ans;

