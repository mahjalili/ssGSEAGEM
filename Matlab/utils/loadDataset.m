% Load one of the datasets.
% A conceptual framework for transcriptional data integration into a genome-scale metabolic model.
%
% INPUTS
%       datasetname - 'GIMME', 'ssGSEA'
%       organism - 'Yeast'
%       dataname - 'Grigaitis', 'GSE8895'
%
% OUTPUTS
%       datasets - reaction expression levels
%
% Author: Mahdi Jalili, 2021

function datasets = loadDataset(datasetname, organism, dataname)
    clear datasets;
    datasets = cell(0);
    data = load(['datasets/',dataname,'/',datasetname,'/expRxns',organism,'.mat']);
    fields = fieldnames(data.expRxns);
    for i = 1:length(fields)
        datasets{i} = struct('name', fields{i}, 'expRxnsMinMax', data.expRxns.(fields{i}).MinMax);
    end
end
