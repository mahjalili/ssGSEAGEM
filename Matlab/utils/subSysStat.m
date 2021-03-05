% Create Subsystem Stat (subSysStatYeast.mat)
% A conceptual framework for transcriptional data integration into a genome-scale metabolic model.
%
% Author: Mahdi Jalili, 2021

clear all;
organism = 'Yeast';

load(['models/',organism,'.mat']);
[subSystems,singleList] = findSubSystemsFromGenes(model);
subSystemsGenes = ({});
subSystemsGenes{length(subSystems),1} = [];
for i = 1:length(singleList)
    for j = 1:length(singleList{i,2})
        idx = find(strcmp(subSystems, singleList{i,2}{j}));
        if idx
            geneidx = find(strcmp(model.genes, singleList{i,1}));
            if ~isempty(model.genes{geneidx,1})
                col = length(find(~cellfun('isempty',subSystemsGenes(idx,:))))+1;
                subSystemsGenes{idx, col} = model.genes{geneidx,1};
            end
        end
    end
end
load('datasets/GIMME/Data.mat');
subSystemsStat = struct();
subSystemsStat.subsystems = subSystems';
subSystemsStat.samples = Data.ColNames;
for i = 1:size(subSystemsGenes,1)
    idx = find(~cellfun('isempty',subSystemsGenes(i,:)));
    genes = unique(subSystemsGenes(i,idx), 'stable');
    subSystemsStat.genes{i,1} = genes;
    [~,ia,~] = intersect(Data.RowNames, genes, 'stable');
    if ia
        for s = 1:length(Data.ColNames)
            subSystemsStat.min{i,s} = min(double(Data(ia, s)));
            subSystemsStat.max{i,s} = max(double(Data(ia, s)));
            subSystemsStat.mean{i,s} = mean(double(Data(ia, s)));
        end
    end
end

save(['datasets/ssGSEA/subSysStat',organism,'.mat'], 'subSystemsStat');
