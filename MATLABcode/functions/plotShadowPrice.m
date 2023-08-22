function plotShadowPrice(model1,model2,labels)
%UNTITLED2 Summary of this function goes here

% identified difference of common mets on shadaw price between model1 and model2

%% Calculate shadow prices for models

OCR=[model1.lb(find(ismember(model1.rxns, 'EX_o2[e]'))):...
    model1.ub(find(ismember(model1.rxns, 'EX_o2[e]')))]';

% model1
model1=changeObjective(model1, 'ATPS4mi');
atpmRates1=zeros(15,1);

for k=1:length(model1.mets);
    ShadowPrices1.('metNames'){k,1}=model1.metNames{k,1};
    % first convert number to string
    k=num2str(k);
    ShadowPrices1.(strcat('met',k))=zeros(15,1);
    % save metabolite name for reference
end

for i=1:length(OCR);
    model1=changeRxnBounds(model1,'EX_o2[e]',OCR(i),'b');
    FBAsolution=optimizeCbModel(model1,'max', 1e-6);
    atpmRates1(i,1)=FBAsolution.f;
    if FBAsolution.stat==1
        % save complete fluxes and shadow prices
        Fluxes{i,1}=FBAsolution.x;
        ShadowPricesTotal{i,1}=FBAsolution.y;
        % save individual shadow prices for metabolites
        for k=1:length(model1.mets)
            % need an index to find the shadow price, otherwise
            % won't work since l is not numerical
            MetInd=find(ismember(model1.mets,model1.mets(k)));
            % first convert number to string
            m=num2str(k);
            % Every value below abs(1e-6) should be considered as 0.
            if abs(FBAsolution.y(MetInd,1))> abs(1e-6);
                ShadowPrices1.(strcat('met',m))(i,1)=FBAsolution.y(MetInd,1);
            else
                ShadowPrices1.(strcat('met',m))(i,1)=0;
            end
        end
    end
end
% All metabolites in model1
allMetabolites1 = model1.mets;
allMetTitles1 = strcat(model1.metNames);

% for model2
model2=changeObjective(model2, 'ATPS4mi');
atpmRates2=zeros(15,1);

for lPD=1:length(model2.mets);
    ShadowPrices2.('metNames'){lPD,1}=model2.metNames{lPD,1};
    % first convert number to string
    lPD=num2str(lPD);
    ShadowPrices2.(strcat('met',lPD))=zeros(15,1);
    % save metabolite name for reference
end

for i=1:length(OCR);
    model2=changeRxnBounds(model2,'EX_o2[e]',OCR(i),'b');
    FBAsolution=optimizeCbModel(model2,'max', 1e-6);
    atpmRates2(i,1)=FBAsolution.f;
    if FBAsolution.stat==1
        % save complete fluxes and shadow prices
        FluxesPD{i,1}=FBAsolution.x;
        ShadowPricesTotalPD{i,1}=FBAsolution.y;
        % save individual shadow prices for metabolites
        for lPD=1:length(model2.mets)
            % need an index to find the shadow price, otherwise
            % won't work since l is not numerical
            MetInd=find(ismember(model2.mets,model2.mets(lPD)));
            % first convert number to string
            n=num2str(lPD);
            % Every value below abs(1e-6) should be considered as 0.
            if abs(FBAsolution.y(MetInd,1))> abs(1e-6);
                ShadowPrices2.(strcat('met',n))(i,1)=FBAsolution.y(MetInd,1);
            else
                ShadowPrices2.(strcat('met',n))(i,1)=0;
            end
        end
    end
end

% All metabolites in model2
allMetabolites2 = model2.mets;
allMetTitles2 = strcat(model2.metNames);


%% Extract common metabolites of interest

% all Common metabolites between control and PD
commonMetabolites = intersect(model1.mets, model2.mets);
commonMetTitles={};
for i = 1:length(commonMetabolites);
    met = find(ismember(model1.mets, commonMetabolites(i)));
    commonMetTitles(end+1,1) = strcat(model1.metNames(met));
end

% Common extracellular metabolites
extracellularMetabolites = {};
exMetTitles = {};
for k=1:length(commonMetabolites)
    confirmEx=findstr(commonMetabolites{k,1}, '[e]');
    if ~isempty(confirmEx)
        m=num2str(k);
        extracellularMetabolites(end+1,1) = commonMetabolites(k);
        exMetTitles(end+1,1)=strcat(commonMetTitles(k));
    end
end

% commom Cytosolic metabolites
cytosolicMetabolites = {};
cytMetTitles = {};
for k=1:length(commonMetabolites)
    confirmEx=findstr(commonMetabolites{k,1}, '[c]');
    if ~isempty(confirmEx)
        m=num2str(k);
        cytosolicMetabolites(end+1,1) = commonMetabolites(k);
        cytMetTitles(end+1,1)=strcat(commonMetTitles(k));
    end
end

metabos.commonMetabolites=commonMetabolites;
metabos.extracellularMetabolites=extracellularMetabolites;
metabos.cytosolicMetabolites=cytosolicMetabolites;

metaboTitles.commonMetTitles=commonMetTitles;
metaboTitles.exMetTitles=exMetTitles;
metaboTitles.cytMetTitles=cytMetTitles;
%% Plot shadow prices for extracellular metabolites

% Insert metabolites of interest (e.g., cytosolic, extracellular, or all)
types=fieldnames(metabos);
titles=fieldnames(metaboTitles);

for i=1:length(types)
    metabolites = metabos.(types{i});
    metaboliteTitles = metaboTitles.(titles{i});
    
    % Find out which metabolite shadow prices are worth plotting. We identify
    % the shadow prices at OCR = -9.53, because this is the inclination point
    % of metabolism in the PD model (SYNPD1).
    allC = find(ismember(model1.mets, metabolites));
    allX = find(ismember(model2.mets, metabolites));
    
    spC = [];
    spX = [];
    for j = 1:length(allC);
        spC(end+1,1) = ShadowPrices1.(strcat('met', num2str(allC(j))))(9,1);
        spX(end+1,1) = ShadowPrices2.(strcat('met', num2str(allX(j))))(9,1);
    end
    summaryTable = table(model1.mets(allC), spC, spX, abs(spC)-abs(spX),abs(abs(spC)-abs(spX)),...
        abs(spC)+abs(spX));
    
    % From a visual analysis of the table, it seems that the top 10 metabolites
    % and the bottom 2 metabolites show the most differences in shadow prices
    % between control and PD. The remaining shadow prices are in the decimal
    % places.
    summaryTable = sortrows(summaryTable,'Var5','descend');
    
    summaryTable.Properties.VariableNames{1} = 'metabolite';
    summaryTable.Properties.VariableNames{4} = 'difference';
    summaryTable.Properties.VariableNames{5} = 'abs_difference';
    
    
    % Delete the unnecessary rows accordingly and extract the new metabolite list:
    summaryTable(13:end,:) = [];
    mets = table2cell(summaryTable(:,1));
    
    b = figure('units','normalized','outerposition',[0 0 1 1]);
    for k=1:length(mets);
        C = find(ismember(model1.mets, mets(k)));
        Cm = num2str(C);
        X = find(ismember(model2.mets, mets(k)));
        Xm=num2str(X);
        data = [ShadowPrices1.(strcat('met', Cm)), ShadowPrices2.(strcat('met', Xm))];
        subplot(3,4,k);
        bar1 = bar(OCR', data);
        set(bar1(1),'FaceColor',[0.2 0.2 1],'BarWidth',1);
        set(bar1(2),'FaceColor',[1 0 0],'BarWidth',1);
        xlabel('Oxygen exchange flux (umol/gDW/hr)','FontSize',12);
        ylabel('Shadow price','FontSize',12);
        t = find(ismember(metabolites, mets(k)));
        title(metaboliteTitles(t),'FontSize',12);
        hold on
    end
    legend({labels{1}, labels{2}}, ...
        'Orientation', 'horizontal',...
        'Position',[0.375 0.0119 0.239 0.0245],...
        'FontSize',12, 'box', 'off');
    % Create textbox
    annotation(b,'textbox',...
        [0.36 0.955 0.27 0.034],...
        'String',{['Shadow prices for ' types{i}]},...
        'FontWeight','bold',...
        'FontSize',14,...
        'EdgeColor','none');
    
    %% Find the minimum and maximum fluxes for a set of reactions of interest
    % robustnessAnalysis_FVA
    
end
end
