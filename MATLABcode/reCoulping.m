function [newmodel] = reCoulping(model)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
modelConstrained=model;
%%% Phosphatidylcholine:
multipleRxnList={'PCHOLP_hs', 'PLA2_2', 'SMS'};
rxnInd = findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1,1];
d=2.025;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,multipleRxnList,c,d,csense);

%%% Adenosine Monophosphate:
multipleRxnList={'AMPDA', 'NTD7'};
rxnInd = findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1];
d=0.2265;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,multipleRxnList,c,d,csense);

%%% Glutamate:
multipleRxnList={'ASPTA', 'ILETA', 'LEUTA', 'VALTA',...
    'ALATA_L', 'GLUCYS', 'GLUDxm', 'GLUDym'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[-1,-1,-1,-1,-1,1,1,1];
d=1.2225;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Aspartate:
multipleRxnList={'ARGSS', 'ASPTA'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1];
d=1.1925;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Serine:
multipleRxnList={'GHMT2r', 'r0060'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1];
d=0.8625;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Arginine:
multipleRxnList={'GLYAMDTRc', 'ARGN', 'r0145'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1,1];
d=0.7245;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Tyrosine:
multipleRxnList={'HMR_6728', 'HMR_6874', 'TYR3MO2', 'TYRTA'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1,1,1];
d=0.55875;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Histidine: 
multipleRxnList={'HISDC', 'HISD'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1];
d=1.095;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Leucine + Isoleucine: 
multipleRxnList={'LEUTA', 'ILETA'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1];
d=1.305;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Valine + Methionine: 
multipleRxnList={'VALTA', 'METAT'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,1];
d=0.705;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

%%% Glycine: 
multipleRxnList={'GTHS', 'GHMT2r'};
rxnInd=findRxnIDs(modelConstrained, multipleRxnList);
c=[1,-1];
d=0.7725;
csense='G';
modelConstrained=constrainRxnListAboveBound(modelConstrained,...
    multipleRxnList,c,d,csense);

newmodel=modelConstrained;
end

