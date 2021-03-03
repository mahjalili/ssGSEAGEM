function contextModel = call_method(model, method, dataset, reserverxns)

%fprintf([dataset.name '(' method ') is in progress ... ']);

fh = str2func([method '.call_' method]);
contextModel = fh(model, dataset, reserverxns);

%fprintf(2, 'Done.\n');

end