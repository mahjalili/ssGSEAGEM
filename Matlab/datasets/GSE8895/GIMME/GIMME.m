[num,txt,~] = xlsread('datasets/GSE8895/GSE8895.xlsx');
nummean = [mean(num(:,1:3),2),mean(num(:,4:6),2),mean(num(:,7:9),2),mean(num(:,10:12),2)];
col = {'Glucose', 'Ethanol', 'Acetate', 'Maltose'};
Data = bioma.data.DataMatrix(nummean, txt(2:end,1), col);
clear num nummean txt col;
%%
% Mean the same genes row
[C, ia, ic] = unique(Data.RowNames, 'stable');
for i=1:length(C)
    idx = find(ismember(ic, i));
    temp = Data(idx(1,1), :);
    temp(1,:) = mean(Data(idx, :), 1);
    if(i==1)
        DataUnique = temp;
    else
        DataUnique = vertcat(DataUnique, temp);
    end
    fprintf(1, [num2str(i) '- ' C{i} ' done.\n']);
end 
Data = DataUnique;
clear C ia ic i idx temp DataUnique;
%%
%Normalize
Data = (Data-min(Data)) ./ (max(Data)-min(Data));
save(['datasets/GSE8895/GIMME/Data.mat'], 'Data');
%%
organism = 'Yeast';
load(['models/',organism,'.mat']);
clear expRxns;
dataset = struct();
dataset.genes = Data.RowNames;
dataset.attr = 'genes';
for i = 1:length(Data.ColNames)
    expField = matlab.lang.makeValidName(Data.ColNames{i});
    fprintf(1, [num2str(i) '- ' expField ' start ...']);
    dataset.genesexp = double(Data(:, i));
    dataset.name = expField;
    expRxns.(expField).MinMax = reactionLevels(model, dataset, @min, @max);
    fprintf(2, ' done.\n');
end
save(['datasets/GSE8895/standard/expRxns',organism,'.mat'], 'expRxns');
%%
dmwrite(Data, 'datasets/GSE8895/GIMME/Data.csv', 'Delimiter', ',');
