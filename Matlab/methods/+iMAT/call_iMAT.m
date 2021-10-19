% Create iMAT Model
%
% INPUTS
%       model - cobra model
%       dataset - 'Grigaitis', 'GSE8895'
%       reserverxns - reserved reactions
%
% OUTPUTS
%       modelMethod - context model
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

function modelMethod = call_iMAT(model, dataset, reserverxns)

    tol = 1e-8;
    runtime = 7200;
    expressionRxns = dataset.expRxnsMinMax;
    threshold_ub = quantile(expressionRxns(expressionRxns~=-1), 0.75);
    threshold_lb = quantile(expressionRxns(expressionRxns~=-1), 0.25);
    expressionRxns(isnan(expressionRxns)) = -1;
    
    modelMethod = iMAT.solve_iMAT(model, expressionRxns, threshold_lb, threshold_ub, tol, reserverxns, '', runtime);

end