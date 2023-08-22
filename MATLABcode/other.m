high_increase=metaPD.total(metaPD.total.inconsistent == 0 & metaPD.total.increased_highfreq,:);
high_increase=sortrows(high_increase,'increased_Realfrequency','descend');

high_decrease=metaPD.total(metaPD.total.inconsistent == 0 & metaPD.total.decreased_highfreq,:);
high_decrease=sortrows(high_decrease,'decreased_Realfrequency','descend');

high_inconsistent=metaPD.total(metaPD.total.inconsistent == 1,:);
high_inconsistent=sortrows(high_inconsistent,'decreased_Realfrequency','descend');

ProteinFile='~/drive/bioenergeticsPD/fromXi/data/Proteomics/PMID: 33276480';
processedProteome=readtable([ProteinFile filesep 'processed_Proteome.xlsx']);
% control
data_control=table2array(processedProteome(:,6:10));

% meanValue = mean(data_control, 2);  % Calculate the mean along each gene (row)
% stdValue = std(data_control, 0, 2); % Calculate the standard deviation along each gene (row)
% normalizedMatrix = (data_control - meanValue) ./ stdValue;

% log2 scale
scaledMatrix_control = log2(data_control);
x=mean(scaledMatrix_control,2)
histogram(x)

controlProtein=[table(processedProteome.EntrezID),table(x)];
controlProtein.Properties.VariableNames={'genes','expVal'};

writetable(controlProtein,[ProteinFile filesep 'controlProtein.csv'])

scaledMean_control=mean(scaledMatrix_control,2);
scaledMean_control=table(scaledMean_control);
scaledMean_control.Properties.RowNames=processedProteome.ProteinIDs;

y1=scaledMean_control(scaledMean_control.scaledMean_control>20,:);

%%
genes=regexprep(model.genes, '\.\d', '');
ID=processedProteome.ProteinIDs(ismember(processedProteome.EntrezID,double(string(genes))));


sum(ismember(y1.Properties.RowNames,ID))

% PD
data_pd=table2array(processedProteome(:,12:16));
scaledMatrix_pd = log2(data_pd);
x2=mean(scaledMatrix_pd,2)
histogram(x2)
title('Log2(LFQ value) PD')
xlabel('mean of each gene')
ylabel('number of genes')

PDprotein=[table(processedProteome.EntrezID),table(x2)];
PDprotein.Properties.VariableNames={'genes','expVal'};

writetable(PDprotein,[ProteinFile filesep 'PDprotein.csv'])


scaledMean_pd=mean(scaledMatrix_pd,2);
scaledMean_pd=table(scaledMean_pd);
scaledMean_pd.Properties.RowNames=processedProteome.ProteinIDs;

y2=scaledMean_pd(scaledMean_pd.scaledMean_pd>12,:);

sum(ismember(y2.Properties.RowNames,ID))

entrez1=processedProteome.EntrezID(ismember(processedProteome.ProteinIDs,y1.Properties.RowNames));
entrez2=processedProteome.EntrezID(ismember(processedProteome.ProteinIDs,y2.Properties.RowNames));

%%% RNA-seq

selected1=SCTnormcontrol(log(SCTnormcontrol.rowMeansNormalize_data_control)>-0.25,:);
selected2=SCTnormpd(log(SCTnormpd.rowMeansNormalize_data_pd)>0,:);
selected3=scaledcontrol(log(scaledcontrol.rowMeansScaled_data_control)>0,:);
selected4=scaledpd(log(scaledpd.rowMeansScaled_data_pd)>0,:);

gene1=selected1.entrze(ismember(selected1.entrze,double(string(genes))))
gene2=selected2.entrze(ismember(selected2.entrze,double(string(genes))))
gene3=selected3.entrze(ismember(selected3.entrze,double(string(genes))))
gene4=selected4.entrze(ismember(selected4.entrze,double(string(genes))))

sum(ismember(entrez1,gene1))
sum(ismember(entrez1,gene2))
sum(ismember(entrez1,gene3))
sum(ismember(entrez1,gene4))

sum(ismember(entrez2,gene1))
sum(ismember(entrez2,gene2))
sum(ismember(entrez2,gene3))
sum(ismember(entrez2,gene4))
%%
load('~/drive/bioenergeticsPD/fromXi/results/fastcoreModels/SYN.mat')
load('~/drive/bioenergeticsPD/fromXi/results/fastcoreModels/SYNPD.mat')
load('~/drive/bioenergeticsPD/fromXi/results/fastcoreModels/ASYN.mat')
load('~/drive/bioenergeticsPD/fromXi/results/fastcoreModels/ASYNPD.mat')

SYN.SYN2=GeneratedModel;
SYN.SYN2Unconstrained=GeneratedModel;
SYNPD.SYNPD2=GeneratedModel;
SYNPD.SYNPD2Unconstrained=GeneratedModel;

ASYN.ASYN2=GeneratedModel;
ASYN.ASYN2Unconstrained=GeneratedModel;
ASYNPD.ASYNPD2=GeneratedModel;
ASYNPD.ASYNPD2Unconstrained=GeneratedModel;
