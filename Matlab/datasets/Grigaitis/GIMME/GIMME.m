[num,txt,~] = xlsread('datasets/Grigaitis/EcoYeast_Transcriptomics_Normalized.xlsx');
nummean = [mean(num(:,1:2),2),mean(num(:,3:4),2),mean(num(:,5:6),2),mean(num(:,7:10),2),mean(num(:,11:13),2),mean(num(:,14:15),2)];
col = {'MA_020','MA_023','MA_027','MA_030','MA_032','MA_034'};
Data = bioma.data.DataMatrix(nummean, txt(2:end,1), col);
%Normalize
Data = (Data-min(Data)) ./ (max(Data)-min(Data));
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
save(['datasets/Grigaitis/GIMME/expRxns',organism,'.mat'], 'expRxns');
%%
dmwrite(Data, 'datasets/Grigaitis/GIMME/Data.csv', 'Delimiter', ',');
