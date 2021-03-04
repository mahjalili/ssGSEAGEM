inits;
organism = 'Yeast';
dataname = 'Grigaitis'; % Grigaitis, GSE8895
%%
%'biomass pseudoreaction': 'r_4041';
%'sce00190  Oxidative phosphorylation': {'r_0099';'r_0226';'r_0227';'r_0438';'r_0439';'r_0531';'r_0568';'r_0569';'r_0721';'r_1021';'r_2142';'r_2143';'r_2144';'r_1085';'r_1086'}
reserverxns = {'r_4041','r_0099','r_0226','r_0227','r_0438','r_0439','r_0531','r_0568','r_0569','r_0721','r_1021','r_2142','r_2143','r_2144','r_1085','r_1086'};
datasets = {'GIMME','ssGSEA'};
methods = {'GIMME'};
load(['models/',organism,'.mat']);
model0 = model; clear model;
for dset = 1:length(datasets)
    data = loadDataset(datasets{dset}, organism, dataname);
    for mthd = 1:length(methods)
        for d = 1:length(data)
            fprintf(1, [datasets{dset} ' ' methods{mthd} ' ' num2str(d) ' start ...']);
            modelp = properModel(model0, dataname, data{d}.name);
            model = call_method(modelp, methods{mthd}, data{d}, reserverxns);
            save(['Results/',dataname,'/',datasets{dset},'/',data{d}.name,'.mat'], 'model');
            fprintf(2, ' done.\n');
        end
    end
end
