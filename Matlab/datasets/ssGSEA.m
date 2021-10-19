%
% Weigth the genes expression with ssGSEA scores and map to model reactions
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

clear all;
load('Matlab/models/Yeast.mat');
datasets = {'Grigaitis', 'GSE8895'};
for ds = 1:length(datasets)
    load(['Matlab/datasets/', datasets{ds}, '/GIMME/Data.mat']);
    score = readtable(['Matlab/ssGSEA/', datasets{ds}, '_ssGSEA-Genescore.csv'], 'ReadRowNames',true);
    for i=1:length(Data.RowNames)
        idx = find(strcmp(score.Properties.RowNames, Data.RowNames{i}));
        if ~isempty(idx)
            if length(idx)>1
                error('more single idx');
            else
                d = double(Data(i,:));
                s = double(score{idx,:});
                Data(i,:) = d .* s;
            end
        end
    end
    DatassGSEA = bioma.data.DataMatrix(Data, Data.RowNames, Data.ColNames);
    save(['Matlab/datasets/', datasets{ds}, '/ssGSEA/DatassGSEA.mat'], 'DatassGSEA');
    %%
    clear expRxns;
    dataset = struct();
    dataset.genes = DatassGSEA.RowNames;
    dataset.attr = 'genes';
    for i = 1:length(DatassGSEA.ColNames)
        expField = DatassGSEA.ColNames{i};
        fprintf(1, [num2str(i) '- ' expField ' start ...']);
        dataset.genesexp = double(DatassGSEA(:, i));
        dataset.name = expField;
        expRxns.(expField).MinMax = reactionLevels(model, dataset, @min, @max);
        fprintf(2, ' done.\n');
    end
    save(['Matlab/datasets/', datasets{ds}, '/ssGSEA/expRxns.mat'], 'expRxns');
end