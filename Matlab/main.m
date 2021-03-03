inits;

type = {'Recon1','Recon2','Recon3D'};    %'Recon1','Recon2','Recon3D'
dataset = 'ssGSEAw';  %'NCI60','ssGSEAw';
methods = {'GIMME'};    % 'GIMME','iMAT','INIT','FASTCORE'
%%
for i = 1:length(type)
    growthEval(dataset, type{i}, methods);
end
%%
for i = 1:length(type)
    findEssGenes(dataset, type{i}, methods);
end
%%
for i = 1:length(type)
    OncoTsNumerator(dataset, type{i}, methods);
end
