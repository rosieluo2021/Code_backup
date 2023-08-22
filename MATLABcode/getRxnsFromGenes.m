function [ReactionList] = getRxnsFromGenes(model, listGenes,MatrixFormatOutput)
%%
% This function gets reactions from a list of genes and return the rxn abr
% This function may be redundant to findRxnsFromGenes.
%
% [ReactionList] = GetRxnsFromGenes(model, listGenes)
% Input
% model                 model structure
% listGenes             cell array of genes, no transcripts (but works on human network as it removes the transcripts id's in model.genes
% MatrixFormatOutput    Optional: output is given in matrix format (default: 1)
% ReactionList      array of rxn abr for each gene in list (first entry corresponds to query gene)
%
% April 2012 IT
% April 2014 IT added option for matrix output

if ~exist('MatrixFormatOutput','var')
    MatrixFormatOutput = 1;
end


G = regexp(model.genes,'\.','split');
for i = 1: length(G)
    G2 = G{i,1};
    Genes(i,1) = G2(1);
end

clear R
T = find(ismember(Genes, listGenes));
if MatrixFormatOutput == 1
    R=cell(length(T),100);
    for i = 1 : length(T)
        X= model.rxns(find(model.rxnGeneMat(:,T(i))));
        R(i,1) = Genes(T(i));
        R(i,2:length(X)+1) = X';
    end
    ReactionList = R;
else
    ReactionList = [];
    for i = 1 : length(T)
        X= model.rxns(find(model.rxnGeneMat(:,T(i))));
        ReactionList = [ReactionList; X];
    end
    ReactionList = unique(ReactionList);
end

