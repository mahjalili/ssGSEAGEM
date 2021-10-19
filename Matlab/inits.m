%
% Initialization
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

clear all; clc;

COBRASOLVER = 'ibm_cplex'; 
changeCobraSolver(COBRASOLVER, 'all');
% Add paths
addpath(genpath('datasets'));
addpath('methods');
addpath(genpath('utils'));

clear COBRASOLVER ans;

