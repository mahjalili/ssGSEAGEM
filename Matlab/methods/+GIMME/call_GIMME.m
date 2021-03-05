% Create GIMME Model
% A conceptual framework for transcriptional data integration into a genome-scale metabolic model.
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

function modelMethod = call_GIMME(model, dataset, reserverxns)

    obj_frac = 0.9;
 
    %For 'A and B'->min(A,B) and 'A or B'->max(A,B)
    expressionRxns = dataset.expRxnsMinMax;
    threshold = quantile(expressionRxns(expressionRxns~=-1), 0.25);
    expressionRxns(isnan(expressionRxns)) = -1;

    [modelMethod, ~] = GIMME.solve_GIMME(model, expressionRxns, threshold, obj_frac, reserverxns);

end
