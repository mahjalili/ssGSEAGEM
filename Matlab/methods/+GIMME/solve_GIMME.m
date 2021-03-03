function [tissueModel, gimmeSolution] = solve_GIMME(model, expressionRxns, threshold, obj_frac, reserverxns)
% Implementation from COBRA toolbox 2017
% Use the GIMME algorithm (`Becker and Palsson, 2008`) to extract a context
% specific model using data. GIMME minimizes usage of low-expression
% reactions while keeping the objective (e.g., biomass) above a certain
% value. Note that this algorithm does not favor the inclusion of reactions
% not related to the objective.
%
% USAGE:
%
%    tissueModel = GIMME(model, expressionRxns, threshold)
%
% INPUTS:
%
%    model:               input model (COBRA model structure)
%    expressionRxns:      expression data, corresponding to model.rxns (see
%                         mapGeneToRxn.m)
%    threshold:           expression threshold, reactions below this are
%                         minimized
%
% OPTIONAL INPUT:
%    obj_frac:            minimum fraction of the objective(s) of model
%                         (default value - 0.9)
%
% OUTPUTS:
%    tissueModel:         extracted model
%
% `Becker and Palsson (2008). Context-specific metabolic networks are
% consistent with experiments. PLoS Comput. Biol. 4, e1000082.`
%
% .. Author: - Originally written by Becker and Palsson, adapted by S. Opdam and A. Richelle - May 2017



    if nargin < 4 || isempty(obj_frac)
        obj_frac =0.9;
    end

    objectiveCol = [find(model.c) obj_frac]; 

    nRxns = size(model.S,2);

    %first make model irreversible
    [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);

    nbExpressionRxns = size(expressionRxns,1);
    if (nbExpressionRxns < nRxns)
        disp('Warning: Fewer expression data inputs than reactions');
        expressionRxns(nbExpressionRxns+1:nRxns,:) = zeros(nRxns-nbExpressionRxns, size(expressionRxns,2));
    end

    nIrrevRxns = size(irrev2rev,1);
    expressionRxnsIrrev = zeros(nIrrevRxns,1);
    for i=1:nIrrevRxns
        expressionRxnsIrrev(i,1) = expressionRxns(irrev2rev(i,1),1);
    end

    nObjectives = size(objectiveCol,1);
    for i=1:nObjectives
        objectiveColIrrev(i,:) = [rev2irrev{objectiveCol(i,1),1}(1,1) objectiveCol(i,2)];
    end

    %Solve initially to get max for each objective
    for i=1:size(objectiveCol)
        %define parameters for initial solution
        modelIrrev.c=zeros(nIrrevRxns,1);
        modelIrrev.c(objectiveColIrrev(i,1),1)=1;

        %find max objective
        FBAsolution = optimizeCbModel(modelIrrev);
        if (FBAsolution.stat ~= 1)
            not_solved=1;
            display('Failed to solve initial FBA problem');
            return
        end
        maxObjective(i)=FBAsolution.f;
    end

    model2gimme = modelIrrev;
    model2gimme.c = zeros(nIrrevRxns,1);


    for i=1:nIrrevRxns
        if (expressionRxnsIrrev(i,1) > -1)   %if not absent reaction
            if (expressionRxnsIrrev(i,1) < threshold)
                model2gimme.c(i,1) = threshold-expressionRxnsIrrev(i,1); %FIX: use expression level as weight
            end
        end
    end

    for i=1:size(objectiveColIrrev,1)
        model2gimme.lb(objectiveColIrrev(i,1),1) = objectiveColIrrev(i,2) * maxObjective(i);
    end

    gimmeSolution = optimizeCbModel(model2gimme,'min');

    if (gimmeSolution.stat ~= 1)
    %No solution for the problem
        display('Failed to solve GIMME problem'); 
        gimmeSolution.x = zeros(nIrrevRxns,1);
    end

    reactionActivityIrrev = zeros(nIrrevRxns,1);
    for i=1:nIrrevRxns
        if ((expressionRxnsIrrev(i,1) > threshold) || (expressionRxnsIrrev(i,1) == -1))
            reactionActivityIrrev(i,1)=1;
        elseif (gimmeSolution.x(i,1) > 0)
            reactionActivityIrrev(i,1)=2;
        end
    end

    %Translate reactionActivity to reversible model
    reactionActivity = zeros(nRxns,1);
    for i=1:nRxns
        for j=1:size(rev2irrev{i,1},2)
            if (reactionActivityIrrev(rev2irrev{i,1}(1,j)) > reactionActivity(i,1))
                reactionActivity(i,1) = reactionActivityIrrev(rev2irrev{i,1}(1,j));
            end
        end
    end
    
    reactionActivity(find(ismember(model.rxns, reserverxns)), 1) = 1;
    remove = model.rxns(reactionActivity == 0);
    removefinal = remove;
    %{
    removefinal = {};
    j = 1;
    for i=1:length(remove)
       dummymodel = removeRxns(model,remove{i});
       sol = optimizeCbModel(dummymodel, 'max');
       if(sol.f > 0)
           removefinal{j} = remove{i};
           j = j+1;
       end
    end
    %}
    tissueModel = removeRxns(model,removefinal); 
end

function [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model, varargin)
% Implementation from COBRA toolbox 2017
% Converts model to irreversible format, either for the entire model or for
% a defined list of reversible reactions.
%
% USAGE:
%    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model, varargin)
%
% INPUT:
%    model:         COBRA model structure
%
% OPTIONAL INPUTS:
%    varargin:      Additional Options as ParameterName/Value pairs.
%                   Allowed Parameters are:
%                   * sRxns: List of specific reversible reactions to convert to
%                     irreversible (Default = model.rxns)
%                   * flipOrientation:  Alter reactions that can only carry negative
%                     flux by flipping and marking them with '_r'
%                     (Default = true)
%                   * orderReactions: Order Reactions such that reverse
%                     reactions directly follow their forward reaction.
% OUTPUTS:
%    modelIrrev:    Model in irreversible format
%    matchRev:      Matching of forward and backward reactions of a reversible reaction
%    rev2irrev:     Matching from reversible to irreversible reactions
%    irrev2rev:     Matching from irreversible to reversible reactions
%
% EXAMPLE:
%
%    % To convert a whole model to an irreversible model:
%    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model)
%
%    % To Convert only Reactions R1, R2 and R15:
%    reactionsToRevert = {'R1', 'R2', 'R15'};
%    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model, 'sRxns', reactionsToRevert)
%
%    % To order the reactions such that backward follows forward reaction:
%    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model, 'orderReactions', true)
%
% NOTE:
%   Uses the reversible list to construct a new model with reversible
%   reactions separated into forward and backward reactions.  Separated
%   reactions are appended with '_f' and '_b' and the reversible list tracks
%   these changes with a '1' corresponding to separated forward reactions.
%   Reactions entirely in the negative direction will be reversed and
%   appended with '_r'.
%
% .. Authors:
%       - written by Gregory Hannum 7/9/05
%       - Modified by Markus Herrgard 7/25/05
%       - Modified by Jan Schellenberger 9/9/09 for speed.
%       - Modified by Diana El Assal & Fatima Monteiro 6/2/17 allow to
%         optionally only split a specific list of reversible reactions to
%         irreversible, without appending '_r'.
%       - Modified by Thomas Pfau June 2017 - Also include all fields
%         associated to reactions and add additional options.


    parser = inputParser();
    parser.addRequired('model',@isstruct);
    parser.addParamValue('sRxns',model.rxns,@iscell);
    parser.addParamValue('flipOrientation',true,@(x) islogical(x) || isnumeric(x));
    parser.addParamValue('orderReactions',false,@(x) islogical(x) || isnumeric(x));

    parser.parse(model,varargin{:});
    sRxns = parser.Results.sRxns;
    flipOrientation = parser.Results.flipOrientation;
    orderReactions = parser.Results.orderReactions;

    %Flip all pure backward reactions and append a _r
    backReacs = ismember(model.rxns,sRxns) & model.lb < 0 & model.ub <= 0;

    if flipOrientation
        model.S(:,backReacs) = -model.S(:,backReacs);
        templbs = -model.ub(backReacs);
        model.ub(backReacs) = -model.lb(backReacs);
        model.lb(backReacs) = templbs;
        model.c(backReacs) = - model.c(backReacs); %Also flip the objective coefficient, as otherwise the target changes.
        model.rxns(backReacs) = strcat(model.rxns(backReacs),'_r');
    end

    %
    % Note: reactions which can only carry negative flux, will have an inactive
    % forward reaction.
    revReacs = ismember(model.rxns,sRxns) & model.lb < 0 & model.ub > 0;
    nRevRxns = sum(revReacs);
    nRxns = numel(model.rxns);
    rxnIDs = 1:nRxns;
    irrevRxnIDs = nRxns + (1:nRevRxns);


    %teat special fields: S, lb, ub, rxns
    model.S(:,end+1:end+nRevRxns) = -model.S(:,revReacs);

    %update the lower and upper bounds (first for the reversed reactions, as
    %otherwise the information is lost
    model.ub(end+1:end+nRevRxns) = max(0,-model.lb(revReacs));
    model.lb(end+1:end+nRevRxns) = max(0,-model.ub(revReacs));
    model.lb(revReacs) = max(0,model.lb(revReacs));
    model.ub(revReacs) = max(0,model.ub(revReacs));

    %Extend the c vector by the negative (otherwise the objective changes)
    model.c(end+1:end+nRevRxns) = -model.c(revReacs);

    %Alter the reaction ids (as defined)
    RelReacNames = model.rxns(revReacs);
    model.rxns(revReacs) = strcat(RelReacNames,'_f');
    model.rxns(end+1:end+nRevRxns) = strcat(RelReacNames,'_b');

    %And update all other relevant fields (which have not yet been altered)
    fields = getModelFieldsForType(model,'rxns','fieldSize',nRxns);
    for i = 1:length(fields)    
        cfield = fields{i};
        if size(model.(cfield),1) == nRxns
            model.(cfield)(end+1:end+nRevRxns,:) = model.(cfield)(revReacs,:);
        elseif size(model.(cfield),2) == nRxns
            model.(cfield)(:,end+1:end+nRevRxns) = model.(cfield)(:,revReacs);
        end
    end


    %Now, map the reactions
    irrev2rev = [rxnIDs';rxnIDs(revReacs)'];
    rev2irrev = num2cell(rxnIDs');
    rev2irrev(revReacs) = num2cell([rxnIDs(revReacs)',irrevRxnIDs'],2);
    matchRev = zeros(size(model.rxns));
    matchRev(revReacs) = irrevRxnIDs;
    matchRev(irrevRxnIDs) = rxnIDs(revReacs);

    reacPosSet = false(numel(model.rxns),1);
    reacorder = zeros(length(model.rxns),1);
    %If a specific order is required (e.g. OptKnock assumes this ordering).
    if orderReactions
        pos = 1;
        newMatchRev = zeros(size(model.rxns));
        for i = 1:size(reacPosSet,1)
            if reacPosSet(i)
                continue;
            end
            reacPosSet(i) = true;
            reacorder(pos) = i;
            irrev2rev(pos) = i;
            rev2irrev{i} = pos;                
            pos = pos + 1;
            if matchRev(i) ~= 0
                %Update the position vector, and the indicator vector.
                newMatchRev(pos-1) = pos;
                newMatchRev(pos) = pos -1;
                irrev2rev(pos) = i;
                rev2irrev{i} = [rev2irrev{i}, pos];
                reacorder(pos) = matchRev(i);
                pos = pos + 1;            
                reacPosSet(matchRev(i)) = true;
            else

            end
        end
        %Now, get the relevant fields for this models rxns again and reorder
        %them....
        fields = getModelFieldsForType(model,'rxns');
        nIrrevRxns = length(model.rxns);
        for i = 1:length(fields)
            cfield = fields{i};
            if size(model.(cfield),1) == nIrrevRxns
                model.(cfield)(:,:) = model.(cfield)(reacorder,:);
            elseif size(model.(cfield),2) == nIrrevRxns
                model.(cfield)(:,:) = model.(cfield)(:,reacorder);
            end
        end
        matchRev = newMatchRev;
    end
    %Mark the model type.
    modelIrrev = model;
    modelIrrev.match = matchRev;
    modelIrrev.reversibleModel = false;    
end