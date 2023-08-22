function [TableChecks, Table_csources, CSourcesTestedRxns, TestSolutionNameOpenSinks,TestSolutionNameClosedSinks] = performSanityChecks_bioenergeticsPD(model,ExtraCellCompIn,ExtraCellCompOut,runSingleGeneDeletion,method)
%from function performSanityChecksonRecon.m
%   made come changes

if ~exist('ExtraCellCompIn','var')
    ExtraCellCompIn = '[e]'; % [e] compartment by default
end
if ~exist('ExtraCellCompOut','var')
    ExtraCellCompOut = '[e]'; % [e] compartment by default
end

if ~exist('runSingleGeneDeletion','var')
    runSingleGeneDeletion = 0; % do not run single gene deletion by default
end

if ~exist('method','var')
    method = 'FBA'; % do not run single gene deletion by default
else
    param.solver = 'mosek';
    param.internalNetFluxBounds='original';
    param.method='fluxes';
    param.printLevel=2;
    param.debug = 1;
    param.feasTol=1e-7;
end


cnt = 1;
tol = 1e-6;

model.rxns(find(ismember(model.rxns,'ATPM')))={'DM_atp_c_'};
model.rxns(find(ismember(model.rxns,'ATPhyd')))={'DM_atp_c_'};
% adds DM_atp to model if not exist

if isempty(strmatch('DM_atp_c_',model.rxns,'exact'))
    [model, rxnIDexists] = addReaction(model,'DM_atp_c_', 'reactionFormula', 'h2o[c] + atp[c]  -> adp[c] + h[c] + pi[c] ');
   % [model, rxnIDexists] = addReaction(model,'DM_atp_c_', 'reactionFormula', '1 atp[c] ->');
end



TestSolutionNameOpenSinks ='';
TestSolutionNameClosedSinks = '';

%% model produces energy from water
if 1
    %%model produces energy from water!
    modelClosed = model;
    
    modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
    modelexchanges4 = strmatch('EX_',modelClosed.rxns);
    modelexchanges2 = strmatch('DM_',modelClosed.rxns);
    modelexchanges3 = strmatch('sink_',modelClosed.rxns);
    selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
    
    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
    modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
    modelClosedATP = changeObjective(modelClosed,'DM_atp_c_');
    modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c_',0,'l');
    modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c_',1000,'u');
    modelClosedATP = changeRxnBounds(modelClosedATP,strcat('EX_h2o',ExtraCellCompIn),-1,'l');
    switch method
        case 'FBA'
            FBA3=optimizeCbModel(modelClosedATP);
        case 'eFBA'
            try
                [FBA3,~]=entropicFluxBalanceAnalysis(modelClosedATP,param);
                FBA3.f=sum(FBA.v(find(model.c>0)));
            catch
                FBA3.f=0;
            end
             if ~isfield(FBA3,'f')
                FBA3.f=0;
            end
            
    end
    
    TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, except h2o';
    if abs(FBA3.f) > 1e-6
        TableChecks{cnt,2} = 'model produces energy from water!';
    else
        TableChecks{cnt,2} = 'model DOES NOT produce energy from water!';
    end
    cnt = cnt + 1;
end
%% model produces energy from water and oxygen!
if 1
    modelClosed = model;
    modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
    modelexchanges4 = strmatch('EX_',modelClosed.rxns);
    modelexchanges2 = strmatch('DM_',modelClosed.rxns);
    modelexchanges3 = strmatch('sink_',modelClosed.rxns);
    selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
    
    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
    modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
    modelClosedATP = changeObjective(modelClosed,'DM_atp_c_');
    modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c_',0,'l');
    modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c_',1000,'u');
    modelClosedATP = changeRxnBounds(modelClosedATP,strcat('EX_h2o',ExtraCellCompIn),-1,'l');
    modelClosedATP = changeRxnBounds(modelClosedATP,strcat('EX_o2',ExtraCellCompIn),-1,'l');
    
    switch method
        case 'FBA'
            FBA6=optimizeCbModel(modelClosedATP);
        case 'eFBA'
            try
                FBA6=entropicFluxBalanceAnalysis(modelClosedATP,param);
                FBA6.f=sum(FBA.v(find(model.c>0)));
            catch
                FBA6.f=0;
            end
             if ~isfield(FBA6,'f')
                FBA6.f=0;
            end
    end
    
    TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, except h2o and o2';
    if abs(FBA6.f) > 1e-6
        TableChecks{cnt,2} = 'model produces energy from water and oxygen!';
    else
        TableChecks{cnt,2} = 'model DOES NOT produce energy from water and oxygen!';
    end
    cnt = cnt + 1;
end
%% model produces matter when atp demand is reversed!
if 1
    modelClosed = model;
    modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
    modelexchanges4 = strmatch('EX_',modelClosed.rxns);
    modelexchanges2 = strmatch('DM_',modelClosed.rxns);
    modelexchanges3 = strmatch('sink_',modelClosed.rxns);
    selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
    
    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
    modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
    
    modelClosed = changeObjective(modelClosed,'DM_atp_c_');
    modelClosed.lb(find(ismember(modelClosed.rxns,'DM_atp_c_'))) = -1000;
    modelClosed.ub(selExc)=1000;
    FBA = optimizeCbModel(modelClosed,'min');
    switch method
        case 'FBA'
            FBA = optimizeCbModel(modelClosed,'min');
        case 'eFBA'
            try
                FBA=entropicFluxBalanceAnalysis(modelClosed,'min',param);
                FBA.f=sum(FBA.v(find(model.c>0)));
            catch
                FBA.f=0;
            end
             if ~isfield(FBA,'f')
                FBA.f=0;
            end
    end
    
    TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, allow DM_atp_c_ to be reversible';
    if abs(FBA.f) > 1e-6
        TableChecks{cnt,2} = 'model produces matter when atp demand is reversed!';
    else
        TableChecks{cnt,2} = 'model DOES NOT produce matter when atp demand is reversed!';
    end
    cnt = cnt + 1;
end
%% model has flux through h[m] demand !
if 1
    modelClosed = model;
    modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
    modelexchanges4 = strmatch('EX_',modelClosed.rxns);
    modelexchanges2 = strmatch('DM_',modelClosed.rxns);
    modelexchanges3 = strmatch('sink_',modelClosed.rxns);
    selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
    
    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
    modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
    modelClosed = addDemandReaction(modelClosed,'h[m]');
    modelClosed = changeObjective(modelClosed,'DM_h[m]');
    modelClosed.ub(find(ismember(modelClosed.rxns,'DM_h[m]'))) = 1000;
    modelClosed.ub(selExc)=1000;
    switch method
        case 'FBA'
            FBA = optimizeCbModel(modelClosed,'max');
        case 'eFBA'
            try
                FBA=entropicFluxBalanceAnalysis(modelClosed,'max');
                FBA.f=sum(FBA.v(find(model.c>0)));
            catch
                FBA.f=0;
            end
             if ~isfield(FBA,'f')
                FBA.f=0;
            end
    end
    TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, test flux through DM_h[m] (max)';
    if abs(FBA.f) > 1e-6
        TableChecks{cnt,2} = 'model has flux through h[m] demand (max)!';
    else
        TableChecks{cnt,2} = 'model has NO flux through h[m] demand (max)!';
    end
    cnt = cnt + 1;
end
%% model has flux through h[c] demand !
if 1
    modelClosed = model;
    modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
    modelexchanges4 = strmatch('EX_',modelClosed.rxns);
    modelexchanges2 = strmatch('DM_',modelClosed.rxns);
    modelexchanges3 = strmatch('sink_',modelClosed.rxns);
    selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
    
    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
    modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
    modelClosed = addDemandReaction(modelClosed,'h[c]');
    modelClosed = changeObjective(modelClosed,'DM_h[c]');
    modelClosed.ub(find(ismember(modelClosed.rxns,'DM_h[c]'))) = 1000;
    modelClosed.ub(selExc)=1000;
    switch method
        case 'FBA'
            FBA = optimizeCbModel(modelClosed,'max');
        case 'eFBA'
            try
            FBA=entropicFluxBalanceAnalysis(modelClosed);
            FBA.f=sum(FBA.v(find(model.c>0)));
            catch
                FBA.f=0;
            end
            if ~isfield(FBA,'f')
                FBA.f=0;
            end
    end
    TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, test flux through DM_h[c] (max)';
    if abs(FBA.f) > 1e-6
        TableChecks{cnt,2} = 'model has flux through h[c] demand (max)!';
    else
        TableChecks{cnt,2} = 'model has NO flux through h[c] demand (max)!';
    end
    cnt = cnt + 1;
end
if 1
    modelClosed = model;
    modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
    modelexchanges4 = strmatch('EX_',modelClosed.rxns);
    modelexchanges2 = strmatch('DM_',modelClosed.rxns);
    modelexchanges3 = strmatch('sink_',modelClosed.rxns);
    selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
    
    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
    modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
    modelClosed = addDemandReaction(modelClosed,'h[c]');
    modelClosed = changeObjective(modelClosed,'DM_h[c]');
    modelClosed.ub(find(ismember(modelClosed.rxns,'DM_h[c]'))) = 1000;
    modelClosed.lb(find(ismember(modelClosed.rxns,'DM_h[c]'))) = -1000;
    modelClosed.ub(selExc)=1000;
    modelClosed.ub(find(ismember(modelClosed.rxns,strcat('EX_h',ExtraCellCompOut)))) = 0;

    switch method
        case 'FBA'
            FBA = optimizeCbModel(modelClosed,'min');
        case 'eFBA'
            try
                FBA=entropicFluxBalanceAnalysis(modelClosed,'min');
                FBA.f=sum(FBA.v(find(model.c>0)));
            catch
                FBA.f=0;
            end
            if ~isfield(FBA,'f')
                FBA.f=0;
            end
    end
    TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, ub of EX_h[e] = 0, test flux through DM_h[c] (min)';
    if abs(FBA.f) > 1e-6
        TableChecks{cnt,2} = 'model has flux through h[c] demand (min)!';
    else
        TableChecks{cnt,2} = 'model has NO flux through h[c] demand (min)!';
    end
    cnt = cnt + 1;
end
%% Test Human Functions
TableChecks{cnt,1} = 'Test metabolic objective functions with open sinks';
if 1 % perform test function
    switch method
        case 'FBA'
            [TestSolution,TestSolutionNameOpenSinks, TestedRxnsSinks, PercSinks] = Test4HumanFctExtv5(model,'all');
            TableChecks{cnt,2} = strcat('Done_FBA. See variable TestSolutionNameOpenSinks for results. The model passes ', num2str(length(find(abs(TestSolution)>tol))),' out of ', num2str(length(TestSolution)), 'tests');
        case 'eFBA'
            param.debug = 1;
            param.feasTol=1e-7;
            optionSinks=0;
            [TestSolution,TestSolutionNameOpenSinks, TestedRxnsSinks, PercSinks] = Test4HumanFctExtv5_eFBA(model,'all',optionSinks,param);
            TableChecks{cnt,2} = strcat('Done_eFBA. See variable TestSolutionNameOpenSinks for results. The model passes ', num2str(length(find(abs(TestSolution)>tol))),' out of ', num2str(length(TestSolution)), 'tests');
    end
else
    TableChecks{cnt,2} = 'Not performed.';
end
cnt =  cnt + 1;

%% Compute ATP yield
TableChecks{cnt,1} = 'Compute ATP yield';
if 1 % test ATP yield
    [Table_csources, CSourcesTestedRxns, Perc] = testATPYieldFromCsources(model,method,ExtraCellCompIn,ExtraCellCompOut);
    TableChecks{cnt,2} = 'Done. See variable Table_csources for results.';
else
    TableChecks{cnt,2} = 'Not performed.';
    CSourcesTestedRxns = '';
end
cnt = cnt + 1;
%%
TableChecks{cnt,1} = 'Check empty columns in rxnGeneMat';
if 1
    E = find(sum(model.rxnGeneMat)==0);
    if isempty(E)
        TableChecks{cnt,2} = 'No empty columns in rxnGeneMat.';
    else
        TableChecks{cnt,2} = 'Empty columns in rxnGeneMat.';
    end
else
    TableChecks{cnt,2} = 'Not performed.';
end
cnt = cnt + 1;
%%
TableChecks{cnt,1} = 'Check that demand reactions have a lb >= 0';
if 1
    DMlb = find(model.lb(strmatch('DM_',model.rxns))<0);
    if isempty(DMlb)
        TableChecks{cnt,2} = 'No demand reaction can have flux in backward direction.';
    else
        TableChecks{cnt,2} = 'Demand reaction can have flux in backward direction.';
    end
else
    TableChecks{cnt,2} = 'Not performed.';
end
cnt = cnt + 1;

%%
TableChecks{cnt,1} = 'Check whether singleGeneDeletion runs smoothly';
if runSingleGeneDeletion == 1
    try
        [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
        TableChecks{cnt,2} = 'singleGeneDeletion finished without problems.';
    catch
        TableChecks{cnt,2} = 'There are problems with singleGeneDeletion.';
    end
else
    TableChecks{cnt,2} = 'Not performed.';
end

end

