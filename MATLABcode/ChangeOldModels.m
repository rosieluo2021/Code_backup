function [newmodel] = ChangeOldModels(model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

modelnames=fieldnames(model);
index=find(contains(lower(string(modelnames)),'1'));
for i=1:length(index)
    modelChange=model.(modelnames{index(i)});
    if isfield(modelChange,'csense') & length(modelChange.csense)~= length(modelChange.mets)
        modelChange.csense((end-10):end)=[];
        newmodelIndex=find(contains(lower(string(modelnames)),'2'));
        Standardmodel=model.(modelnames{newmodelIndex(1)});
        modelChange.dsense=Standardmodel.dsense;
        modelChange=findSExRxnInd(modelChange);
    end
    if ~isfield(modelChange,'C') & ~contains(modelnames{index(i)},'Unconstrained')
        try
        modelChange=reCoulping(modelChange);
        catch ME
         modelChange.ctrs=Standardmodel.ctrs;
         modelChange.d=Standardmodel.d;
         modelChange.dsense=Standardmodel.dsense;
        end
        modelChange=findSExRxnInd(modelChange);
    end
    if isfield(modelChange,'C') & islogical(modelChange.C)
        modelChange.C=sparse(double(modelChange.C));
    end
    model.(modelnames{index(i)})=modelChange;
end

newmodel=model;
end

