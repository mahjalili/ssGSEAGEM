%
% Proper Model
%
% INPUTS
%       model - cobra model
%       dataname - 'Grigaitis', 'GSE8895'
%       cond - database condition
%
% OUTPUTS
%       pmodel - propered model
%
% Author: Mahdi Jalili, 2021
% Mahdi Jalili, Pranas Grigaitis, Martin Scharm, Olaf Wolkenhauer, and Ali Salehzadeh-Yazdi. Metabolic function-based normalization improves transcriptome data-driven reduction of genome-scale metabolic models.

function [pmodel] = properModel(model, dataname, cond)
    %Biomass
    model.c(find(model.c)) = 0; % growth
    model.c(3413) = 1;  % biomass pseudoreaction
    %lb
    if dataname == "Grigaitis"
        MA_020 = [2.18762135908281,5.48243868451356,0,5.56477430729108];
        MA_023 = [2.53191609379646,5.96729153529064,0,5.95237174521335];
        MA_027 = [3.17989324245154,7.29368953491503,0,7.86045670199146];
        MA_030 = [4.25742952315959,6.79030452659990,2.24412644654585,7.97114466109714];
        MA_032 = [4.54631700656014,6.64257299155141,3.40100677035195,8.39496225260453];
        MA_034 = [5.82503479656592,5.00481124983181,7.90649585743139,8.06639755874342];
        rxns = [2588,2815,2629,2549];
    elseif dataname == "GSE8895"
        Glucose	= [2.74,2.85,1.15,0,0,0];
        Maltose	= [3.05,3.05,0,0.61,0,0];
        Ethanol	= [6.87,3.26,0,0,3.78,0];
        Acetate	= [7.4,7.45,0,0,0,5.89];
        rxns = [2815,2549,2588,2776,2629,2519];
        rxnsO = [2174,2445,2447,2448,2453,2454,2533,2564,2580,2583,2584,2589,2592];
        for i=1:length(rxnsO)
            model.lb(rxnsO(i)) = 0;
        end
    end
    lb = eval(cond);
    for i=1:length(rxns)
        model.lb(rxns(i)) = lb(i) * -1;
    end
    pmodel = model;
end
