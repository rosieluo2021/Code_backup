%clear all
clear
beep on
initCobraToolbox

%% 1. genericModel
%genericModelName = 'Recon3DModel_301.mat'; % recon3DModel
genericModelName = 'Recon3d.mat' ;% Preconditional recon model

switch genericModelName
    case 'Recon3DModel_301.mat'
        load('~/fork-cobratoolbox/test/models/mat/Recon3DModel_301.mat');
        model=Recon3D;
    case 'Recon3d.mat'
        %load('~/drive/bioenergeticsPD/rawData/Recon3d.mat')
        load('~/drive/bioenergeticsPD/fromXi/Recon3d.mat')
        modelOrig=model;
end
%% 2. Context-specific data
dataFolder = ['~/drive/bioenergeticsPD/fromXi/reogranizeData'];

%%%% add GEO and nonGEO results in the reconstruction file
% bibliomicData = 'synaptic_bibliomicData.xlsx';
% bibliomicData = 'synapticPD_bibliomicData.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData_unconstrained.xlsx';
% bibliomicData = 'astro_bibliomicData.xlsx';
% bibliomicData = 'astro_bibliomicData_unconstrained.xlsx';
% specificData = preprocessingOmicsModel([dataFolder filesep bibliomicData], 1, 1);

%%%% reconstruction file remove corresponding inactive genes with nonGEO
% bibliomicData = 'synaptic_bibliomicData_New.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_New.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained_New.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained_New.xlsx';
%specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep bibliomicData], 1, 1);

%%%% reconstruction file remove corresponding inactive genes (without GEO and nonGEO
% bibliomicData = 'synaptic_bibliomicData_New1.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_New1.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained_New1.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained_New1.xlsx';
% specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new1' filesep bibliomicData], 1, 1);

%%%% raw reconstruction file from Diana (with nonGEO active genes)
% bibliomicData = 'synaptic_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_unconstrained_reconstruction.xlsx';
bibliomicData = 'nonsynapticPD_bibliomicData_unconstrained_reconstruction.xlsx';
specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new2' filesep bibliomicData], 1, 1);

%%%% only raw reconstruction file from Diana (without GEO and nonGEO active genes)
% bibliomicData = 'synaptic_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData_reconstruction.xlsx';
% bibliomicData = 'nonsynaptic_bibliomicData_unconstrained_reconstruction.xlsx';
% bibliomicData = 'nonsynapticPD_bibliomicData_unconstrained_reconstruction.xlsx';
% specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new3' filesep bibliomicData], 1, 1);

%%%% raw reconstruction file from Diana (without GEO and nonGEO active
%%%% genes) + newly added active rxns
% bibliomicData = 'synaptic_bibliomicData_reconstruction1.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_reconstruction1.xlsx';
% bibliomicData = 'synaptic_bibliomicData_unconstrained_reconstruction1.xlsx';
% bibliomicData = 'synapticPD_bibliomicData_unconstrained_reconstruction1.xlsx';

% specificData = preprocessingOmicsModel([dataFolder filesep 'new' filesep 'new4' filesep bibliomicData], 1, 1);
%% use GEO transcriptomic expression value
dataFile='~/drive/bioenergeticsPD/fromXi/data/RNA-seq/GSE157783/data';
if any(strfind(bibliomicData,'synapticPD')) & ~any(strfind(bibliomicData,'astro'))
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allPD_geneval.txt']);
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allPD_zscore.txt']);
%     specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);
% snRNA-seq
    specificData.transcriptomicData = readtable([dataFile filesep 'SCTnorm_pd.csv']);
    %specificData.transcriptomicData = readtable([dataFile filesep 'scaled_pd.csv']);
    param.transcriptomicThreshold=-0.5;
elseif ~any(strfind(bibliomicData,'synapticPD')) & ~any(strfind(bibliomicData,'astro'))
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allControl_geneval.txt']);
%     specificData.transcriptomicData = readtable([dataFolder filesep 'GEO_transcriptomicData_allControl_zscore.txt']);
%     specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);
% snRNA-seq
    specificData.transcriptomicData = readtable([dataFile filesep 'SCTnorm_control.csv']);
    %specificData.transcriptomicData = readtable([dataFile filesep 'scaled_control.csv']);
    param.transcriptomicThreshold=-0.5;
end
%% use proteomic data from Sarah Plum et al. 2020
ProteinFile='~/drive/bioenergeticsPD/fromXi/data/Proteomics/PMID: 33276480';

if any(strfind(bibliomicData,'PD')) & ~any(strfind(bibliomicData,'astro'))
    specificData.proteomicData=readtable([ProteinFile filesep 'PDprotein.csv']);
    param.thresholdP=22;%22 20
    
elseif  ~any(strfind(bibliomicData,'PD')) & ~any(strfind(bibliomicData,'astro'))
    specificData.proteomicData=readtable([ProteinFile filesep 'controlProtein.csv']);
    param.thresholdP=21;%21 18
end
%% 3.Technical parameters
%choose solver
solver='ibm_cplex';
solverOK = changeCobraSolver(solver,'LP');
% parameters
param.TolMinBoundary = -1e3;
param.TolMaxBoundary =  1e3;
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 100; % fluxEpsilon=1e-4;
param.closeIons = false; %not specified
param.closeUptakes = false; %not specified
% param.nonCoreSinksDemands = 'closeNone';
param.nonCoreSinksDemands = 'closeAll';
param.activeGenesApproach = 'oneRxnPerActiveGene';
% param.tissueSpecificSolver = 'fastCore'; % Choose fastcore
param.tissueSpecificSolver = 'thermoKernel'; % Choose thermoKernel
param.fluxEpsilon = feasTol * 100;% fluxEpsilon=1e-4;
param.fluxCCmethod = 'fastcc';
param.addCoupledRxns = 0; % have been added in the preconditional recon model
param.thermoFluxEpsilon=param.fluxEpsilon;
param.weightsFromOmics=true;

param.inactiveGenesTranscriptomics = 1; 

%% 4.XomicsToModel pipeline to generate models
modelname=regexprep(bibliomicData,['_(\w+)ta'],'');
modelname=regexprep(modelname,'.xlsx','');
modelname=[modelname '_model'];

switch param.tissueSpecificSolver
    case 'fastCore'
resultsFolder = ['~/drive/bioenergeticsPD/fromXi/model/fastcore/new' filesep modelname];
    case 'thermoKernel'
resultsFolder = ['~/drive/bioenergeticsPD/fromXi/model/thermokernel/new' filesep modelname];
end

if ~isfolder(resultsFolder) 
    mkdir(resultsFolder);
end
cd(resultsFolder);

param.printLevel = 1;
param.debug = true;
if isunix()
    name = getenv('USER');
else
    name = getenv('username');
end
param.diaryFilename = [resultsFolder filesep datestr(now,30) '_' name '_diary.txt'];
%%
[GeneratedModel, modelGenerationReport] = XomicsToModel(model, specificData, param);
GeneratedModel.lb(contains(GeneratedModel.rxns,'DM_clpn_hs[c]'))=0.000649;
save([resultsFolder filesep [modelname '.mat']],'GeneratedModel')