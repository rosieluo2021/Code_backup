function [CharacteristicTable] = ModelCharacteristic(models,printFlag)
%To generate a table with the Characteristic results for each model
% USAGE:
%
%
% INPUT
%
% OUTPUT

if ~exist('printFlag', 'var') || isempty(printFlag)
    printFlag = 0;
elseif (~isnumeric(printFlag) & ~islogical(printFlag))
    error('printFlag should be a number or a bool')
end

Newmodels=models;

if isstruct(Newmodels)
    % model characteristics
    Characteristic=cell(10,length(fieldnames(Newmodels)));
    modelname=fieldnames(Newmodels);
    for i=1:length(modelname)
        model=Newmodels.(modelname{i});
        Characteristic{1,i}=length(unique(model.genes));
        Mrxns=findRxnFromCompartment(model,'[m]');
        Characteristic{2,i}=length(findGenesFromRxns(model,Mrxns(:,1)));
        Characteristic{3,i}=length(unique(model.rxns));
        Characteristic{4,i}=size(Mrxns,1);
        Trans={'Transport, mitochondrial','Transport, extracellular','Transport, endoplasmic reticular'...
            'Transport, golgi apparatus','Transport, lysosomal','Transport, nuclear','Transport, peroxisomal'};
        TransRxns=[];
        for j=1:length(Trans)
            rxns=findRxnsFromSubSystem(model,Trans{j});
            TransRxns=[TransRxns ; rxns];
        end
        Characteristic{5,i}=length(unique(TransRxns));
        if ~isfield(model,'ExchRxnBool') || ~isfield(model,'DMRxnBool') || ~isfield(model,'SinkRxnBool')
            model=findSExRxnInd(model);
        end
        Characteristic{6,i}=length(unique(model.rxns(model.ExchRxnBool==1)));
        FBAsolution=optimizeCbModel(model);
        Characteristic{7,i}=length(unique(model.rxns(FBAsolution.v<0 & model.ExchRxnBool)));
        compartments=regexp(model.mets, '\[(.*?)\]', 'match');
        Characteristic{8,i}=length(unique(string(compartments)));
        Characteristic{9,i}=length(unique(regexprep(model.mets,'(\[\w\])','')));
        Characteristic{10,i}=sum(string(compartments)=='[m]');
    end

Characteristic=cell2table(Characteristic);
Characteristic.Properties.VariableNames=fieldnames(Newmodels);
Characteristic.Properties.RowNames={'Total Num of genes','Genes in [m]','Total Num of rxns', 'Mitochondrial rxns','Transport Rxns','Exchange rxns','Uptake rxns (Ex rxns with FBA.v<0)','Compartments','Total Num of unique Mets', 'Metabolites in [m]'};
CharacteristicTable=Characteristic;
else
    disp('please check the input model')
end
%% Print tables with output if printFlag = 1
if printFlag ==1
    CharacteristicTable
end
