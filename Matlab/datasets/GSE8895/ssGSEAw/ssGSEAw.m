clear all;
type = 'Yeast';
%%
data = readtable(['ssGSEA/GSE8895/ssGSEA',type,'.csv'],'ReadRowNames',true);
DatassGSEAw = bioma.data.DataMatrix(data{:,:}, data.Properties.RowNames, data.Properties.VariableNames);
load('datasets/GSE8895/standard/Data.mat');
load(['datasets/GSE8895/ssGSEAw/subSysStat',type,'.mat']);
load(['models/',type,'.mat']);
[~,singleList] = findSubSystemsFromGenes(model);
mat = zeros(length(DatassGSEAw.RowNames),length(DatassGSEAw.ColNames));
for i=1:length(DatassGSEAw.RowNames)
    %idx = find(strcmp(model.genes, DatassGSEAw.RowNames{i}));
    idx = find(strcmp(Data.RowNames, model.genes{i}));
    if isempty(idx)
        mat0 = zeros(1,length(DatassGSEAw.ColNames));mat0(mat0==0)=NaN;
        for j = 1:length(singleList{i,2})
            idx2 = find(strcmp(subSystemsStat.subsystems, singleList{i,2}{j}));
            if idx2 && ~isempty(subSystemsStat.min{idx2})
                disp([num2str(i),'mahdi']);
                mat0(end+1,:) = subSystemsStat.min{idx2,:};
            else
                disp([num2str(i),'mahdi000']);
            end
        end
        mat(i,:) = min(mat0);
    elseif length(idx)>1 
        error('Error');
    else
        mat(i,:) = double(Data(idx,:));
    end
end
    
DatassGSEAw = bioma.data.DataMatrix(DatassGSEAw.*mat, DatassGSEAw.RowNames, DatassGSEAw.ColNames);
save(['datasets/GSE8895/ssGSEAw/DatassGSEAw',type,'.mat'], 'DatassGSEAw');

%%
clear expRxns;
dataset = struct();
dataset.genes = DatassGSEAw.RowNames;
dataset.attr = 'genes';
for i = 1:length(DatassGSEAw.ColNames)
    expField = DatassGSEAw.ColNames{i};
    fprintf(1, [num2str(i) '- ' expField ' start ...']);
    dataset.genesexp = double(DatassGSEAw(:, i));
    dataset.name = expField;
    expRxns.(expField).MinMax = reactionLevels(model, dataset, @min, @max);
    expRxns.(expField).Q25 = quantile(dataset.genesexp, 0.25);
    expRxns.(expField).Q75 = quantile(dataset.genesexp, 0.75);
    fprintf(2, ' done.\n');
end
save(['datasets/GSE8895/ssGSEAw/expRxns',type,'.mat'], 'expRxns');

