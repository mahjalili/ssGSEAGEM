%
% Make Models
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

run('Matlab/inits.m');
%%
datasets = {'Grigaitis', 'GSE8895'};
methods = {'standard', 'ssGSEA'};
integrations = {'GIMME', 'iMAT', 'INIT'};
% Reserve the biomass pseudoreaction (r_4041) and Oxidative phosphorylation
reserverxns = {'r_4041','r_0099','r_0226','r_0227','r_0438','r_0439','r_0531','r_0568','r_0569','r_0721','r_1021','r_2142','r_2143','r_2144','r_1085','r_1086'};
load('Matlab/models/Yeast.mat');
model0 = model; clear model;
%% create the models
for dset = 1:length(datasets)
    for meth = 1:length(methods)
        data = loadDataset(methods{meth}, datasets{dset});
        for integ = 1:length(integrations)
            for d = 1:length(data)
                fprintf(1, [datasets{dset}, '-', methods{meth}, ' ' integrations{integ}, ' ', num2str(d), ' start ...']);
                modelp = properModel(model0, datasets{dset}, data{d}.name);
                model = call_method(modelp, integrations{integ}, data{d}, reserverxns);
                save(['Matlab/Results/',datasets{dset},'/',integrations{integ},'/',methods{meth},'/',data{d}.name,'.mat'], 'model');
                fprintf(2, ' done.\n');
            end
        end
    end
end