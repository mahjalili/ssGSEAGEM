%
% Load one of the datasets.
%
% INPUTS
%       datasetname - 'GIMME', 'ssGSEA'
%       dataname - 'Grigaitis', 'GSE8895'
%
% OUTPUTS
%       datasets - reaction expression levels
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

function datasets = loadDataset(datasetname, dataname)
    clear datasets;
    datasets = cell(0);
    data = load(['Matlab/datasets/',dataname,'/',datasetname,'/expRxns.mat']);
    fields = fieldnames(data.expRxns);
    for i = 1:length(fields)
        datasets{i} = struct('name', fields{i}, 'expRxnsMinMax', data.expRxns.(fields{i}).MinMax);
    end
end
