% Create Contex Model
% A conceptual framework for transcriptional data integration into a genome-scale metabolic model.
%
% INPUTS
%       model - cobra model
%       method - 'GIMME'
%       dataset - 'Grigaitis', 'GSE8895'
%       reserverxns - reserved reactions
%
% OUTPUTS
%       contextModel - context model
%
% Author: Mahdi Jalili, 2021

function contextModel = call_method(model, method, dataset, reserverxns)
   fh = str2func([method '.call_' method]);
  contextModel = fh(model, dataset, reserverxns);
end
