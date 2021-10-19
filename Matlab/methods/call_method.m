%
% INPUTS
%       model - cobra model
%       method - 'GIMME', 'iMAT', 'INIT'
%       dataset - 'Grigaitis', 'GSE8895'
%       reserverxns - reserved reactions
%
% OUTPUTS
%       contextModel - context model
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

function contextModel = call_method(model, method, dataset, reserverxns)
   fh = str2func([method '.call_' method]);
  contextModel = fh(model, dataset, reserverxns);
end
