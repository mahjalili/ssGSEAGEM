%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

function modelrev = addRevField(model)

    % Add depreciated rev field to model
    if ~isfield('rev', model)
        sRxns = model.rxns;
        revReacs = ismember(model.rxns, sRxns) & model.lb < 0 & model.ub > 0;
        model.rev = revReacs;
    end
    modelrev = model;

end