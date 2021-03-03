function contextModel = funcModel( contextModel )

    paramConsistency.epsilon=1e-10;
    paramConsistency.modeFlag=0;
    paramConsistency.method='fastcc';
    
    [~, fluxConsistentRxnBool] = findFluxConsistentSubset(contextModel, paramConsistency);
    remove = contextModel.rxns(fluxConsistentRxnBool==0);
    contextModel = removeRxns(contextModel, remove);
    contextModel = removeUnusedGenes(contextModel);
end