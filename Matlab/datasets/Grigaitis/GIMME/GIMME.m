[num,txt,~] = xlsread('datasets/Grigaitis/EcoYeast_Transcriptomics_Normalized.xlsx');
nummean = [mean(num(:,1:2),2),mean(num(:,3:4),2),mean(num(:,5:6),2),mean(num(:,7:10),2),mean(num(:,11:13),2),mean(num(:,14:15),2)];
col = {'MA_020','MA_023','MA_027','MA_030','MA_032','MA_034'};
Data = bioma.data.DataMatrix(nummean, txt(2:end,1), col);
clear num nummean txt col;
%%
% Mean the same genes row
%no dublicated
%{
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
clear C ia ic i idx temp;
%}
%%
%Normalize
%{
    dim = 2, then zscore uses the means and standard deviations along the rows of X.
    A = magic(4);
    data = zscore(A,0,1);
    normA = data - min(data); normA = normA ./ max(normA);
    %OR
    normA = (A-min(A)) ./ (max(A)-min(A));
%}
Data = (Data-min(Data)) ./ (max(Data)-min(Data));
%%
type = 'Yeast';
%%
%load(['models/',type,'-consistent.mat']);
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
    expRxns.(expField).Q25 = quantile(dataset.genesexp, 0.25);
    expRxns.(expField).Q75 = quantile(dataset.genesexp, 0.75);
    fprintf(2, ' done.\n');
end
save(['datasets/Grigaitis/standard/expRxns',type,'.mat'], 'expRxns');

%%
%Data2 = colnames(Data, [1:6], fieldnames(expRxns));
dmwrite(Data, 'datasets/Grigaitis/standard/Data.csv', 'Delimiter', ',');
%%