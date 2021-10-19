% Create INIT Model
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

function modelMethod = call_INIT(model, dataset, reserverxns)

    tol = 1e-8;
    expressionRxns = dataset.expRxnsMinMax;
    threshold = quantile(expressionRxns(expressionRxns~=-1), 0.25);
    expressionRxns(expressionRxns<=threshold) = -expressionRxns(expressionRxns<=threshold);

    modelMethod = INIT.solve_INIT(model, expressionRxns, tol, reserverxns);

end