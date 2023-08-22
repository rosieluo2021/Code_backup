function  FluxInSubsystems(model1, model2,objective1,objective2)
% This code calculate FBA For output of xomicstomodel and plot the flux of
% all reactions in different pathways
%% changing solver to gurobi
changeCobraSolver('gurobi','all');
format long
%% change objective function of the models to biomass_reaction
model1 = changeObjective(model1,objective1);
model2 = changeObjective(model2,objective2);
%% Running FBA
FBAsolution1 = optimizeCbModel(model1,'max');
FBAsolution2 = optimizeCbModel(model2,'max');
subSystems = intersect(getModelSubSystems(model1),getModelSubSystems(model2));

%% Ploting figures for common subsystems
disp('Ploting figures for common subsystems')
for i= 1:length(subSystems)
    
    SubSystem = subSystems(i);
    rxnID1 = findRxnIDs(model1,findRxnsFromSubSystem(model1,SubSystem));
    rxnID2 = findRxnIDs(model2,findRxnsFromSubSystem(model2,SubSystem));
    %     rxnName = intersect(findRxnsFromSubSystem(model1,SubSystem),findRxnsFromSubSystem(model2,SubSystem));   
    % remove the rxns with flux <=1e-5
    remove1 = [];
    for i = 1:length(rxnID1)
        if (abs(FBAsolution1.v(rxnID1(i))) <= 1e-5)
            FBAsolution1.v(rxnID1(i)) = 0;
            if FBAsolution1.v(rxnID1(i)) == 0
                remove1 = [remove1, rxnID1(i)];
            end
        end
    end
    rxnID1_N = setdiff(rxnID1, remove1);
    remove2 = [];
    for i = 1:length(rxnID2)
        if (abs(FBAsolution2.v(rxnID2(i))) <= 1e-5)
            FBAsolution2.v(rxnID2(i)) = 0;
            if FBAsolution2.v(rxnID2(i)) == 0
                remove2 = [remove2, rxnID2(i)];
            end
        end
    end
    rxnID2_N = setdiff(rxnID2, remove2);
    
    if  isempty(rxnID1_N) & isempty(rxnID2_N)
        fprintf('All fluxes in %s subsystem in model1 and model2 are zero.%s\n', SubSystem{:})
        fprintf('\n')
    elseif ~isempty(rxnID1_N) & isempty(rxnID2_N)
        
        if (FBAsolution1.v(rxnID1_N)) == 0
            fprintf('No fluxes in %s subsystem in model1 are zero.%s\n', SubSystem{:})
            fprintf('\n')
        elseif any((FBAsolution1.v(rxnID1_N)) ~= 0)
            %%%%%%%%%       plot model1
            c = categorical(model1.rxnNames(rxnID1_N));
            c = string(c);
            c = strrep(c,'_','-');
            [m, ~] = size(rxnID1_N);
            b = FBAsolution1.v(rxnID1_N);
            figure
            bar (b,'FaceColor' , [0 0.4470 0.7410])
            set(gca,'XTick',[1:m],'xticklabel',c)
            set(gca,'XTickLabelRotation',45, 'FontSize',12)
            ylabel('Flux(umol/gDW)', 'FontSize',14,'FontWeight','bold')
            title([char(SubSystem) ' in model1'],'FontSize',12);
            legend('model1');
            h = gca;
            h.YAxis.FontWeight = 'bold';
            h.XAxis.FontWeight = 'bold';
        end
    elseif isempty(rxnID1_N) & ~isempty(rxnID2_N)
        if(FBAsolution2.v(rxnID2_N)) == 0
            fprintf('No fluxes in %s subsystem in model2 are zero.%s\n', SubSystem{:})
            fprintf('\n')
        elseif any((FBAsolution2.v(rxnID2_N)) ~= 0)
            %%%%%%%       plot model2
            c = categorical(model2.rxnNames(rxnID2_N));
            c = string(c);
            c = strrep(c,'_','-');
            [m, ~] = size(rxnID2_N);
            figure
            d = FBAsolution2.v(rxnID2_N);
            bar (d, 'FaceColor',[0.8500 0.3250 0.0980])
            set(gca,'XTick',[1:m],'xticklabel',c)
            set(gca,'XTickLabelRotation',45, 'FontSize',12)
            ylabel('Flux(umol/gDW)', 'FontSize',14,'FontWeight','bold')
            title([char(SubSystem) ' in model2'],'FontSize',12);
            legend('model1','model2');
            h = gca;
            h.YAxis.FontWeight = 'bold';
            h.XAxis.FontWeight = 'bold';
        end
    elseif ~isempty(rxnID1_N) & ~isempty(rxnID2_N)
        if any((FBAsolution1.v(rxnID1_N)) ~= 0) & any((FBAsolution2.v(rxnID2_N)) ~= 0)
            %%%%%%%%%%% plot comparison
            rxnID = unique([model1.rxns(rxnID1_N);model2.rxns(rxnID2_N)]);
            a = rxnID;
            for i=1:length(a)
                if ismember(a{i},model1.rxns)
                    c(i) = model1.rxnNames(findRxnIDs(model1,a{i}));%intersect((model1.rxnNames(rxnID1_N)),(model2.rxnNames(rxnID2_N)));
                    c(i) = strrep(c(i),'_','-');
                else
                    c(i) = model2.rxnNames(findRxnIDs(model2,a{i}));%intersect((model1.rxnNames(rxnID1_N)),(model2.rxnNames(rxnID2_N)));
                    c(i) = strrep(c(i),'_','-');
                end
            end
            [m, ~] = size(rxnID);
            y = zeros(m,2);
            for i=1:size(a,1)
                model1ID=findRxnIDs(model1,a(i));
                if model1ID==0
                    y(i,1)=0;
                else
                    y(i,1)= FBAsolution1.v(model1ID);
                end
                model2ID=findRxnIDs(model2,a(i));
                if model2ID==0
                    y(i,2)=0;
                else
                    y(i,2)= FBAsolution2.v(findRxnIDs(model2,a(i)));
                end
            end
            
            if ~any(y~=0)
                fprintf('All of the fluxes in the reactions in %s subsystem are zero.%s\n', SubSystem{:})
            elseif isempty(y)
                fprintf('flux empty')
            else
                figure
                bar(y)
                set(gca,'XTick',[1:m],'xticklabel',c)
                set(gca,'XTickLabelRotation',20, 'FontSize',8)
                ylabel('Flux(umol/gDW)', 'FontSize',10,'FontWeight','bold')
                title([char(SubSystem) ' in two models'],'FontSize',8);
                legend('model1','model2');
                legend('Location','best')
                h = gca;
                h.YAxis.FontWeight = 'bold';
                h.XAxis.FontWeight = 'bold';
            end
        end
    end
end
%% plot flux in subsystems that are not common
disp('Ploting figures for uncommon subsystems')
p_subsystems = setxor(subSystems, getModelSubSystems(model1));
for i= 1:length(p_subsystems)
    SubSystem = p_subsystems(i)
    rxnName = model1.rxns(ismember(model1.subSystems,SubSystem));
    rxnID1_N = findRxnIDs(model1,rxnName);
    if (FBAsolution1.v(rxnID1_N))== 0 
        fprintf('All of the Fluxes in %s in model1 are zero.%s\n', SubSystem{:})
        fprintf('\n')
    else
%%%%%%%%%plot model1
        c = categorical(model1.rxnNames(rxnID1_N));
        c = string(c);
        c = strrep(c,'_','-');
        [m, ~] = size(rxnID1);
        b = FBAsolution1.v(rxnID1_N);
        figure
        bar (b, 'FaceColor' ,[0 0.4470 0.7410])
        set(gca,'XTick',[1:m],'xticklabel',c)
        set(gca,'XTickLabelRotation',45, 'FontSize',12)
        ylabel('Flux(umol/gDW)', 'FontSize',14,'FontWeight','bold')
        title(char(SubSystem),'FontSize',12);
        legend('model1');
        h = gca;
        h.YAxis.FontWeight = 'bold';
        h.XAxis.FontWeight = 'bold';
    end
end
%%%%%%%%%plot Tmodel2

q_subsystems = setxor(subSystems, getModelSubSystems(model2));

for i= 1:length(q_subsystems)
    SubSystem = q_subsystems(i);
    rxnName = model2.rxns(ismember(model2.subSystems,SubSystem));
    rxnID2_N = findRxnIDs(model2,rxnName);
    if (FBAsolution2.v(rxnID2_N)) == 0
        fprintf('All of the Fluxes in %s in model2 are zero.%s\n', SubSystem{:})
        fprintf('\n')
    else
        c = categorical(model1.rxns(rxnID2_N));
        c = string(c);
        c = strrep(c,'_','-');
        [m, ~] = size(rxnID2_N);
        d = FBAsolution1.v(rxnID2_N);
        figure
        bar (d, 'FaceColor',[0.8500 0.3250 0.0980])
        c = string(c);
        set(gca,'XTick',[1:m],'xticklabel',c) 
        set(gca,'XTickLabelRotation',45, 'FontSize',12)
        ylabel('Flux(umol/gDW)', 'FontSize',14,'FontWeight','bold')
        title(char(SubSystem),'FontSize',12);
        legend('model2');
        h = gca;
        h.YAxis.FontWeight = 'bold';
        h.XAxis.FontWeight = 'bold';
    end
end
