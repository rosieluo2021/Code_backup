function plotGlycolysisRobustness(Allmodels,ignoreUnconstrain,method)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~isstruct(Allmodels)
    error('check input models')
end

if ~exist('ignoreUnconstrain','var')
    ignoreUnconstrain=0;
end

if ~exist('method','var')
    method='FBA';
else
    param.debug = 1;
    param.feasTol=1e-7;
end

if any(contains(fieldnames(Allmodels),{'SYNPD'}))
    names = [fieldnames(Allmodels.SYN); fieldnames(Allmodels.SYNPD)];
    Allmodels.allSYN = cell2struct([struct2cell(Allmodels.SYN); struct2cell(Allmodels.SYNPD)], names, 1);
    names = [fieldnames(Allmodels.ASYN); fieldnames(Allmodels.ASYNPD)];
    Allmodels.allASYN = cell2struct([struct2cell(Allmodels.ASYN); struct2cell(Allmodels.ASYNPD)], names, 1);
    
    Allmodels=rmfield(Allmodels,{'SYN','SYNPD','ASYN','ASYNPD'});
end

types=fieldnames(Allmodels);
for i=1:length(types)
    models=Allmodels.(types{i});
    modelnum=fieldnames(models);
    if ignoreUnconstrain ~= 0
        modelnum=modelnum(~contains(modelnum,'constrain'));% ignore unconstrained models
    end
    if ~any(contains(modelnum,'ASTRO'))
        if any(contains(modelnum,'PD'))
            % PD and control
            group.index1=contains(modelnum,'1');
            group.index2=contains(modelnum,'2');
            d = figure('units','normalized','outerposition',[0 0 1 1]);
            for k=1:length(fieldnames(group))
                name=fieldnames(group);
                sets=modelnum(group.(name{k}));
                model1=Allmodels.(types{i}).(sets{1});
                model2=Allmodels.(types{i}).(sets{2});
                OCR = [model1.lb(find(ismember(model1.rxns, 'EX_glc_D[e]'))):0]';
                for m =1:length(OCR);
                    model1 = changeRxnBounds(model1,'EX_glc_D[e]',OCR(m,1),'b');
                    model1 = changeObjective(model1,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                            %CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                            ATPS4mi_ox(m,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                            %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                                %CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                            catch
                                %YOOm3i_ox(k,1)=0;
                                ATPS4mi_ox(m,1)=0;
                            end
                    end
                end
                subplot(1,2,k);
                plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                hold on
                clear m n ATPS4mi_ox
                %%%%%%%%%%%%
                % add another model
                for n =1:length(OCR);
                    model2 = changeRxnBounds(model2,'EX_glc_D[e]',OCR(n,1),'b');
                    model2 = changeObjective(model2,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                            %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                            ATPS4mi_ox(n,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                            %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                                %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                            catch
                                %                             CYOOm3i_ox(k,1)=0;
                                ATPS4mi_ox(n,1)=0;
                            end
                    end
                end
                plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                xlabel('Glucose exchange flux (\mumol/gDW/h)', 'FontSize', 24);
                ylabel('ATP synthase flux (\mumol/gDW/h)','FontSize', 24);
                title([sets{1} ' vs ' sets{2}])
                legend(sets{1},sets{2},'box','off','FontSize', 14)
                clear m n OCR ATPS4mi_ox
            end
            clear group d name sets model1 model2
            %
            %
            % old and new
            group.index1=~contains(modelnum,'PD');
            group.index2=contains(modelnum,'PD');
            d = figure('units','normalized','outerposition',[0 0 1 1]);
            for k=1:length(fieldnames(group))
                name=fieldnames(group);
                sets=modelnum(group.(name{k}));
                model1=Allmodels.(types{i}).(sets{1});
                model2=Allmodels.(types{i}).(sets{2});
                OCR = [model1.lb(find(ismember(model1.rxns, 'EX_glc_D[e]'))):0]';
                for m =1:length(OCR);
                    model1 = changeRxnBounds(model1,'EX_glc_D[e]',OCR(m,1),'b');
                    model1 = changeObjective(model1,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                            %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                            ATPS4mi_ox(m,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                            %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                                %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(m,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                            catch
                                %                             CYOOm3i_ox(k,1)=0;
                                ATPS4mi_ox(m,1)=0;
                            end
                    end
                end
                subplot(1,2,k);
                plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                hold on
                %%%%%%%%%%%%
                % add another model
                for n =1:length(OCR);
                    model2 = changeRxnBounds(model2,'EX_glc_D[e]',OCR(n,1),'b');
                    model2 = changeObjective(model2,'ATPM');
                    switch method
                        case 'FBA'
                            FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                            %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                            ATPS4mi_ox(n,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                            %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                        case 'eFBA'
                            try
                                FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                                %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                                ATPS4mi_ox(n,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                            catch
                                %                             CYOOm3i_ox(k,1)=0;
                                ATPS4mi_ox(n,1)=0;
                            end
                            
                    end
                end
                plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
                xlabel('Glucose exchange flux (\mumol/gDW/h)', 'FontSize', 24);
                ylabel('ATP synthase flux (\mumol/gDW/h)','FontSize', 24);
                title([sets{1} ' vs ' sets{2}])
                legend(sets{1},sets{2},'box','off','FontSize', 14)
                clear m n OCR ATPS4mi_ox
            end
        end
    else
        if ~any(contains(modelnum,'constrain'))
            %old model
            model1 = Allmodels.(types{i}).(modelnum{1});
            % Define the flux ranges based on the lower and upper bounds of the
            % corresponding models
            % Using same bounds for control and PD models.
            OCR = [model1.lb(find(ismember(model1.rxns, 'EX_glc_D[e]'))):0]';
            for k =1:length(OCR);
                model1 = changeRxnBounds(model1,'EX_glc_D[e]',OCR(k,1),'b');
                model1 = changeObjective(model1,'ATPM');
                switch method
                    case 'FBA'
                        FBAenergetics = optimizeCbModel(model1,'max', 1e-6);
                        %                         CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                        ATPS4mi_ox(k,1) = FBAenergetics.x(find(contains(model1.rxns, 'ATPS4mi')));
                        %ATPM_ox(i,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                    case 'eFBA'
                        try
                            FBAenergetics=entropicFluxBalanceAnalysis(model1,param);
                            %                             CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                            ATPS4mi_ox(k,1) = FBAenergetics.v(find(contains(model1.rxns, 'ATPS4mi')));
                        catch
                            %                             CYOOm3i_ox(k,1)=0;
                            ATPS4mi_ox(k,1)=0;
                        end
                end
            end
            d = figure('units','normalized','outerposition',[0 0 1 1]);
            plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
            hold on
            %
            %
            %new model
            model2 = Allmodels.(types{i}).(modelnum{2});
            % Define the flux ranges based on the lower and upper bounds of the
            % corresponding models
            % Using same bounds for control and PD models.
            clear ATPS4mi_ox
            for k =1:length(OCR);
                model2 = changeRxnBounds(model2,'EX_glc_D[e]',OCR(k,1),'b');
                model2 = changeObjective(model2,'ATPM');
                switch method
                    case 'FBA'
                        FBAenergetics = optimizeCbModel(model2,'max', 1e-6);
                        if FBAenergetics.stat ==1
                        %CYOOm3i_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'CYOOm3i')));
                        ATPS4mi_ox(k,1) = FBAenergetics.x(find(contains(model2.rxns, 'ATPS4mi')));
                        %ATPM_ox(k,1) = FBAenergetics.x(find(ismember(model1.rxns, 'ATPM')));
                        else
                            ATPS4mi_ox(k,1) =0;
                        end
                    case 'eFBA'
                        try
                            FBAenergetics=entropicFluxBalanceAnalysis(model2,param);
                            %CYOOm3i_ox(k,1) = FBAenergetics.v(find(ismember(model1.rxns, 'CYOOm3i')));
                            ATPS4mi_ox(k,1) = FBAenergetics.v(find(contains(model2.rxns, 'ATPS4mi')));
                        catch
                            %CYOOm3i_ox(k,1)=0;
                            ATPS4mi_ox(k,1)=0;
                        end
                end
            end
            plot(OCR, ATPS4mi_ox, '-*', 'LineWidth', 2)
            xlabel('Glucose exchange flux (\mumol/gDW/h)', 'FontSize', 24);
            ylabel('ATP synthase flux (\mumol/gDW/h)','FontSize', 24);
            hold off
            title([modelnum{1} ' vs ' modelnum{2}])
            legend(modelnum{1},modelnum{2},'box','off','FontSize', 14)
%             ylim(gca,[2, 7])
            clear ATPS4mi_ox
        end
    end
end



