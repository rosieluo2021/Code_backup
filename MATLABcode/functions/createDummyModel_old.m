function [modelDummy,CoreOrganRxns,WeightsPerRxn] = createDummyModel(model,EntrezGeneID , EntrezGeneIDWeights, TolMaxBoundary)
% I generalized dummy reactions for all genes - IT 14.10.01 -  tested
% added option for AddDummy

%initiate variables
NotInOrganWeight = 50;
NoAssignmentWeights = 10;
ToDelWeight = 200; % to avoid leaky models
ToCorrWeight = 50; % to avoid leaky models
WeightsPerRxn.rxns = '';
WeightsPerRxn.weights = [];
RxnsConsidered = '';
%CoreOrganRxns='';
%proteinMap threshold
threshold = 0.2;

if ~exist('TolMaxBoundary','var')
    TolMaxBoundary=1000;
end

if exist('EntrezGeneIDWeights','var') && ~isempty(EntrezGeneIDWeights)
    % get all proteins NOT detected in organ
    EntrezGeneIDZeroWeights = EntrezGeneID (EntrezGeneIDWeights==0);
    % get corresponding reactions
    modelTmp =model;
    modelTmp.genes = regexprep(modelTmp.genes,'\.(\d+)','');
    genes = modelTmp.genes;
    %   modelTmp.genes = strcat('b',genes);
    %  Ptmp = strcat('b',EntrezGeneIDZeroWeights);
    [NotInOrganRxns] = getRxnsFromGenes(modelTmp,EntrezGeneIDZeroWeights,0);
    % set all these to 0 -- this is too stron
    %model.lb(ismember(model.rxns,NotInOrganRxns))=0;
    %model.ub(ismember(model.rxns,NotInOrganRxns))=0;
    % set reaction weights to NotInOrganWeight
    WeightsPerRxn.rxns = [WeightsPerRxn.rxns;NotInOrganRxns];
    weights = NotInOrganWeight*ones(length(NotInOrganRxns),1);
    WeightsPerRxn.weights = [WeightsPerRxn.weights;weights];
    RxnsConsidered = [RxnsConsidered;NotInOrganRxns];
    % get all proteins detected in organ with more than 0.9 - these will be
    % core proteins and result in core reactions
    CoreOrganProtein = EntrezGeneID (EntrezGeneIDWeights>=threshold);
    % get corresponding reactions
    %  [CoreOrganRxnsMatrix] = getRxnsFromGenes(model, CoreOrganProtein,1);
end

% create a dummy reaction for all genes in model
G = regexp(model.genes,'\.','split');
for i = 1: length(G)
    G2 = G{i,1};
    genes(i,1) = G2(1);
end
[RxnsMatrix] = getRxnsFromGenes(model, genes,1);

% create a dummy reaction for these
cntR = 1;
for j = 1 : size(RxnsMatrix,1)
    Rxns = RxnsMatrix(j,2:end);
    for k=1:length(Rxns)
        if ~isempty(char(Rxns{k}))
            Rxns2{k} = char(Rxns(k));
        end
    end
    if exist('Rxns2','var')
        Rxns = Rxns2;
        clear Rxns2
        R = find(ismember(model.rxns,Rxns));
        % add a dummy row
        % to deal with alternate transcripts we will combine them into 1
        % dummy rxn
        %  for k = 1 :length(R)
        %     dummyMet = strcat('dummy_Rxn_',CoreOrganRxnsMatrix(j,1),'_',Rxns(k));
        dummyMet = strcat('dummy_Rxn_',RxnsMatrix(j,1));
        dummyR = strcat('dummy_Rxns_',RxnsMatrix(j,1));
        if ~isempty(find(ismember(model.mets,dummyMet)))
            row = find(ismember(model.mets,dummyMet));
        else
            row = length(model.mets)+1;
        end
        if ~isempty(find(ismember(model.rxns,dummyR)))
            col = find(ismember(model.rxns,dummyR));
        else
            col = length(model.rxns)+1;
        end
        if isfield(model, 'mets')
            model.mets(row) = dummyMet;
        end
        if isfield(model, 'S')
            model.S(row,R) = 1;
        end
        if isfield(model, 'S')
            model.S(row,col) = -1;
        end
        if isfield(model, 'b')
            model.b(row)=0;
        end
        if isfield(model,'SIntMetBool')
            model.SIntMetBool(row) = 1;
        end
        if isfield(model,'SIntRxnBool')
            model.SIntRxnBool(col) = 0;
        end
        if isfield(model,'SConsistentMetBool')
            model.SConsistentMetBool(row) = 1;
        end
        if isfield(model,'SConsistentRxnBool')
            model.SConsistentRxnBool(col) = 0;
        end
        if isfield(model,'fluxConsistentMetBool')
            model.fluxConsistentMetBool(row) = 1;
        end
        if isfield(model,'fluxConsistentRxnBool')
            model.fluxConsistentRxnBool(col) = 1;
        end
        if isfield(model,'thermoFluxConsistentMetBool')
            model.thermoFluxConsistentMetBool(row) = 1;
        end
        if isfield(model,'thermoFluxConsistentRxnBool')
            model.thermoFluxConsistentRxnBool(col) = 1;
        end
%         if isfield(model, 'metNames')
%             model.metNames(row)=dummyMet;
%         end
%         if isfield(model, 'metFormulas') 
%             model.metFormulas(row)={'N'};
%         end
%         if isfield(model, 'metCharge')
%             model.metCharge(row)=0;
%         end
%         if isfield(model, 'metCHEBIID')
%             model.metCHEBIID(row)={'N'};
%         end
%         if isfield(model, 'metInchiString')
%             model.metInchiString(row)={'M'};
%         end
%         if isfield(model, 'metKeggID')
%             model.metKeggID(row)={'M'};
%         end
%         if isfield(model, 'metPubChemID')
%             model.metPubChemID(row)={'M'};
%         end
%         if isfield(model, 'metCharges')
%             model.metCharges(row)=0;
%         end
%         if isfield(model, 'metSmiles')
%             model.metSmiles(row)={'M'};
%         end
%         if isfield(model, 'metHMDBID')
%             model.metHMDBID(row)={'M'};
%         end
%         if isfield(model, 'metInChIString')
%             model.metInChIString(row)={'M'};
%         end
%         if isfield(model, 'metKEGGID')
%             model.metKEGGID(row)={'M'};
%         end
%         if isfield(model, 'metPdMap')
%             model.metPdMap(row) = {''};
%         end
        if isfield(model, 'rxns')
            model.rxns(col) = dummyR;
        end
        if isfield(model, 'lb')
            model.lb(col)=-TolMaxBoundary;
        end
        if isfield(model, 'ub')
            model.ub(col)=TolMaxBoundary;
        end
%         if isfield(model, 'rev')
%             model.rev(col)=0; 
%         end
        if isfield(model, 'c')
            model.c(col)=0;
        end
        if isfield(model, 'rxnGeneMat')
            model.rxnGeneMat(col,:)=0; 
        end
        if isfield(model, 'rules')
            model.rules{col}='';
        end
%         if isfield(model, 'grRules')
%             model.grRules(col)={''};
%         end
%         % confidence score change to admit double instead of cell
%         % (j.modamio 15.12.2017)model.rxnConfidenceScores(col)={'1'}
%         if isfield(model, 'subSystems')
%             model.subSystems(col)={'Dummy_Rxn'};
%         end
        if isfield(model, 'rxnNames')
            model.rxnNames{col}=model.rxns{col};
        end
        if isfield(model,'subSystems')
            model.subSystems{col}= 'Dummy';
        end
%         if isfield(model, 'rxnNotes')
%             model.rxnNotes(col)={''}; 
%         end
%         if isfield(model, 'rxnReferences') 
%             model.rxnReferences(col)={''};
%         end
%         if isfield(model, 'rxnECNumbers')
%             model.rxnECNumbers(col)={''};
%         end
%         if isfield(model, 'rxnConfidenceScores')
%             model.rxnConfidenceScores(col)='1';
%         end
%         if isfield(model, 'rxnKeggID')
%             model.rxnKeggID(col)={''};
%         end
%         if isfield(model, 'rxnConfidenceEcoIDA')
%             model.rxnConfidenceEcoIDA(col)={''};
%         end
%         if isfield(model, 'rxnsboTerm')
%             model.rxnsboTerm(col)={''};
%         end
%         if isfield(model, 'rxnKEGGID')
%             model.rxnKEGGID(col)={''};
%         end
        % 20190228 Aga: csense and C were not updated, causes verifyModel() to fail in fastcc()
        if isfield(model, 'csense')
            model.csense(row) = 'E';
        end
        if isfield(model, 'ctrs')
            model.C(length(model.ctrs),col) = sparse(1,1);
        end
%         if isfield(model, 'rxnCOG')
%             model.rxnCOG(col) = {''};
%         end
%         if isfield(model, 'rxnKeggOrthology')
%             model.rxnKeggOrthology(col) = {''};
%         end
%         if isfield(model, 'rxnReconMap')
%             model.rxnReconMap(col) = {''};
%         end
%         if isfield(model, 'constraintDescription')
%             model.constraintDescription(col) = {''};
%         end
        % if the gene has been found in CoreOrganProtein then the dummy reaction
        % will be added to core reaction
        % this dummy reaction will be a core reactions
        if exist('EntrezGeneIDWeights','var')
            if ~isempty(find(ismember(CoreOrganProtein,genes(j))))
                CoreOrganRxns(cntR) = model.rxns(col);
                cntR = cntR+1;
            end
        end
        RxnsConsidered = [RxnsConsidered; Rxns'];
    end
end
modelInput=model;


CoreRxnAbbr = CoreOrganRxns;

% assign weights to all reactions that had changed bounds in
% reconRecon
% to ToDelWeight

%    WeightsPerRxn.rxns = [WeightsPerRxn.rxns;ChangeLbUbR];
%    weights = ToDelWeight*ones(length(ChangeLbUbR),1); % more likely than those that have not been found
%    WeightsPerRxn.weights = [WeightsPerRxn.weights;weights];
% corrected direction
% ToCorrWeight
%   WeightsPerRxn.rxns = [WeightsPerRxn.rxns;ChangeLbR;ChangeUbR];
%   weights = ToCorrWeight*ones(length(ChangeLbR)+length(ChangeUbR),1); % more likely than those that have not been found
%   WeightsPerRxn.weights = [WeightsPerRxn.weights;weights];

%   RxnsConsidered = [RxnsConsidered;ChangeUbR;ChangeLbR;ChangeLbUbR];
% assign weights to all reactions not considered by organ information
% to NoAssignmentWeights
NoAssignment = setdiff(model.rxns,RxnsConsidered);
WeightsPerRxn.rxns = [WeightsPerRxn.rxns;NoAssignment];
weights = NoAssignmentWeights*ones(length(NoAssignment),1); % more likely than those that have not been found
WeightsPerRxn.weights = [WeightsPerRxn.weights;weights];

WeightsPerRxn = WeightsPerRxn;

%TODO debug why this is so
%addition of dummy reactions can cause flux inconsistency (e.g.
%{'dummy_Rxns_471'})
paramConsistency.epsilon = getCobraSolverParams('LP', 'feasTol');
paramConsistency.method = 'fastcc';
[~, fluxConsistentRxnBool] = findFluxConsistentSubset(model, paramConsistency);
if any(~fluxConsistentRxnBool)
    model.rxns(fluxConsistentRxnBool);
    warning([int2str(nnz(~fluxConsistentRxnBool)) ' flux inconsistent reaction(s) after dummy model creation, removed.'])
    model.S = model.S(:,fluxConsistentRxnBool);
    model.rxns = model.rxns(fluxConsistentRxnBool);
    model.rxnNames = model.rxnNames(fluxConsistentRxnBool);
    model.lb = model.lb(fluxConsistentRxnBool);
    model.ub = model.ub(fluxConsistentRxnBool);
    model.c = model.c(fluxConsistentRxnBool);
    if isfield(model, 'ctrs')
        model.C = model.C(:,fluxConsistentRxnBool);
    end
    if isfield(model, 'rxnGeneMat')
        model.rxnGeneMat = model.rxnGeneMat(fluxConsistentRxnBool,:);
    end
    if isfield(model, 'rules')
        model.rules = model.rules(fluxConsistentRxnBool);
    end
    if isfield(model,'subSystems')
        model.subSystems = model.subSystems(fluxConsistentRxnBool);
    end
    if isfield(model,'SIntRxnBool')
        model.SIntRxnBool = model.SIntRxnBool(fluxConsistentRxnBool);
    end
    if isfield(model,'SConsistentRxnBool')
        model.SConsistentRxnBool = model.SConsistentRxnBool(fluxConsistentRxnBool);
    end
    if isfield(model,'fluxConsistentRxnBool')
        model.fluxConsistentRxnBool = model.fluxConsistentRxnBool(fluxConsistentRxnBool);
    end
    if isfield(model,'thermoFluxConsistentRxnBool')
        model.thermoFluxConsistentRxnBool = model.thermoFluxConsistentRxnBool(fluxConsistentRxnBool);
    end
end

%any zero rows or columns are considered inconsistent
zeroRowBool=~any(model.S,2);
zeroColBool=~any(model.S,1)';
if any(zeroRowBool) || any(zeroColBool)
     warning('%6u\t%6u\t%s\n',nnz(zeroRowBool),nnz(zeroColBool),' zero rows and columns in model.S of dummy model.')
end

modelDummy = model;

