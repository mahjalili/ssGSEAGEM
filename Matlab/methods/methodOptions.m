% Methods options
% GIMME
    % expressionRxnsMAX       reaction expression, expression data corresponding to model.rxns. Note : If no gene-expression data are available for the reactions, set the value to -1
    % threshold            expression threshold, reactions below this are minimized
    %options.GIMME.threshold = 12;
    % obj_frac            minimum fraction of the model objective function
    options.GIMME.obj_frac = 0.9; % OPTIONAL
% /GIMME
% iMAT
    % expressionRxnsMAX       reaction expression, expression data corresponding to model.rxns. Note : If no gene-expression data are available for the reactions, set the value to -1
    % threshold_lb         lower bound of expression threshold, reactions with expression below this value are "non-expressed"
    options.iMAT.threshold_lb = 10;
    % threshold_ub         upper bound of expression threshold, reactions with expression above this value are "expressed" 
    options.iMAT.threshold_ub = 100;
    % tol                 minimum flux threshold for "expressed" reactions (default 1e-8)
    options.iMAT.tol = 1e-8; % OPTIONAL
    % core                cell with reaction names (strings) that are manually put in the high confidence set (default - no core reactions)
    options.iMAT.core = {}; % OPTIONAL
    % runtime             maximum solve time for the MILP (default - 7200s)
    options.iMAT.runtime = 7200; % OPTIONAL
    % osense
    options.iMAT.osense = 'MAX';
% /iMAT
% INIT
    options.INIT = {};
% /INIT
% FASTCORE
    % epsilon   smallest flux value that is considered nonzero (default 1e-4)
    options.FASTCORE.epsilon = 1e-4; % OPTIONAL
    options.FASTCORE.osense = 'MAX';
% /FASTCORE