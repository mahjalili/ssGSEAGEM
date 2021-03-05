% Create Reaction Expression Levels (expRxnsYeast.mat)
% A conceptual framework for transcriptional data integration into a genome-scale metabolic model.
%
% Author: Mahdi Jalili, 2021

clear all;
organism = 'Yeast';
data = readtable(['ssGSEA/ssGSEA',organism,'.csv'],'ReadRowNames',true);
DatassGSEAw = bioma.data.DataMatrix(data{:,:}, data.Properties.RowNames, data.Properties.VariableNames);
load('datasets/GIMME/Data.mat');
load(['datasets/ssGSEA/subSysStat',organism,'.mat']);
load(['models/',organism,'.mat']);
[~,singleList] = findSubSystemsFromGenes(model);
mat = zeros(length(DatassGSEA.RowNames),length(DatassGSEA.ColNames));
for i=1:length(DatassGSEA.RowNames)
    idx = find(strcmp(Data.RowNames, model.genes{i}));
    if isempty(idx)
        mat0 = zeros(1,length(DatassGSEA.ColNames));mat0(mat0==0)=NaN;
        for j = 1:length(singleList{i,2})
            idx2 = find(strcmp(subSystemsStat.subsystems, singleList{i,2}{j}));
            if idx2 && ~isempty(subSystemsStat.min{idx2})
                mat0(end+1,:) = subSystemsStat.min{idx2,:};
            end
        end
        mat(i,:) = min(mat0);
    elseif length(idx)>1 
        error('Error');
    else
        mat(i,:) = double(Data(idx,:));
    end
end
DatassGSEA = bioma.data.DataMatrix(DatassGSEA.*mat, DatassGSEA.RowNames, DatassGSEA.ColNames);
save(['datasets/ssGSEA/DatassGSEA',organism,'.mat'], 'DatassGSEA');
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
save(['datasets/ssGSEA/expRxns',organism,'.mat'], 'expRxns');
