function datasets = loadDataset(datasetname, organism, dataname)
% Load one of the datasets.
%
    clear datasets;
    datasets = cell(0);
    data = load(['datasets/',dataname,'/',datasetname,'/expRxns',organism,'.mat']);
    fields = fieldnames(data.expRxns);
    for i = 1:length(fields)
        datasets{i} = struct('name', fields{i}, 'expRxnsMinMax', data.expRxns.(fields{i}).MinMax, 'Q25', data.expRxns.(fields{i}).Q25, 'Q75', data.expRxns.(fields{i}).Q75);
    end
end