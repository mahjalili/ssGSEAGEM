%
% Create genes score according ssGSEA result
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

clear all;
dataset = {{'Grigaitis', 6}, {'GSE8895', 4}};
load('Matlab/ssGSEA/geneInGeneset.mat');
for ds = 1:length(dataset)
    ssGSEAresult = importdata(['R/out/', dataset{ds}{1}, '_ssGSEA/', dataset{ds}{1}, '_ssGSEA-combined.csv'], ',', 1);
    genesets = ssGSEAresult.textdata(2:end, 1);
    ssGSEAdata = ssGSEAresult.data;
    samples = dataset{ds}{2};
    out = [];
    for s = 1:samples
        normData = rescale(ssGSEAdata(:,s), 1, 10);
        levels = zeros(length(geneInGeneset),1); levels(levels==0) = NaN;
        for i = 1:length(geneInGeneset)
            for j = 1:length(geneInGeneset{i,2})
                idx = find(strcmp(genesets, geneInGeneset{i,2}{j}));
                if idx
                    if isempty(levels(i))
                        levels(i) = double(normData(idx));
                    else
                        levels(i) = min(levels(i), double(normData(idx)));
                    end
                end
            end
        end
        out{s} = levels;
    end
    conditions = replace(replace(split(ssGSEAresult.textdata(1, 1),','), '"', ''), 'score.', '');
    T = table(out{1:samples}, 'RowNames',geneInGeneset(:,1), 'VariableNames',conditions(1:samples));
    writetable(T, ['Matlab/ssGSEA/', dataset{ds}{1}, '_ssGSEA-GeneScore.csv'], 'WriteRowNames',true);
end