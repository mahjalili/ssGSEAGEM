%
% Load the gene expression data and map to model reactions.
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

datasets = {'Grigaitis', 'GSE8895'};
for ds = 1:length(datasets)
    %% Load data
    if datasets{ds} == "Grigaitis"
        [num,txt,~] = xlsread('Matlab/datasets/Grigaitis/EcoYeast_Transcriptomics_Normalized.xlsx'); %See ref.
        nummean = [mean(num(:,1:2),2),mean(num(:,3:4),2),mean(num(:,5:6),2),mean(num(:,7:10),2),mean(num(:,11:13),2),mean(num(:,14:15),2)];
        col = {'MA_020','MA_023','MA_027','MA_030','MA_032','MA_034'};
        Data = bioma.data.DataMatrix(nummean, txt(2:end,1), col);
    else
        % GSE8895
        [num,txt,~] = xlsread('Matlab/datasets/GSE8895/GSE8895.xlsx');
        nummean = [mean(num(:,1:3),2),mean(num(:,4:6),2),mean(num(:,7:9),2),mean(num(:,10:12),2)];
        col = {'Glucose', 'Ethanol', 'Acetate', 'Maltose'};
        Data = bioma.data.DataMatrix(nummean, txt(2:end,1), col);
        clear num nummean txt col;
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
    end
    %%
    % Normalize
    Data = (Data-min(Data)) ./ (max(Data)-min(Data));
    dmwrite(Data, 'Matlab/datasets/', datasets{ds}, '/GIMME/Data.csv', 'Delimiter', ',');
    %%
    % Map data to reactions
    load('Matlab/models/Yeast.mat');
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
    save(['Matlab/datasets/', datasets{ds}, '/GIMME/expRxns.mat'], 'expRxns');
end