function plotRobustnessResults_OldandNew(demandATP)
%plot old and new
%   change demandATP structure
if any(contains(fieldnames(demandATP),{'SYNPD'}))
    names = [fieldnames(demandATP.SYN); fieldnames(demandATP.SYNPD)];
    demandATP.allSYN = cell2struct([struct2cell(demandATP.SYN); struct2cell(demandATP.SYNPD)], names, 1);
    names = [fieldnames(demandATP.ASYN); fieldnames(demandATP.ASYNPD)];
    demandATP.allASYN = cell2struct([struct2cell(demandATP.ASYN); struct2cell(demandATP.ASYNPD)], names, 1);
    
    demandATP=rmfield(demandATP,{'SYN','SYNPD','ASYN','ASYNPD'});
end
types=fieldnames(demandATP);

for i=1:length(types)
    models=fieldnames(demandATP.(types{i}));
    models=models(~contains(models,'constrain'));% ignore unconstrained models
    % 1. plot for two models comparison
    a = figure('units','normalized','outerposition',[0 0 1 1]);
    if contains(models{1},'ASTRO')
        d = 1:10:50; %reduced to the maximum ATPM in astrocytes
    else
        d=10:50:700;
    end
    
    rxns = {'ATPS4mi', 'glycolysis', 'PGK', 'PYK', 'r0408', ...
        'RE0124C', 'SUCOASm', 'NDPK6', 'NDPK2', 'NMNATr'};
    w1 = 0.6;
    w2 = 0.3;
    if any(contains(models,'PD'))
        % PD and control
        if any(contains(models,'1'))
            type.index1=contains(models,'1');
        end
        if any(contains(models,'2'))
            type.index2=contains(models,'2');
        end
        
        % compare old and new models
        group.index1=contains(models,'PD');
        group.index2=~contains(models,'PD');
        for k=1:length(fieldnames(group))
            name=fieldnames(group);
            sets=models(group.(name{k}));
            model1=demandATP.(types{i}).(sets{1});
            model2=demandATP.(types{i}).(sets{2});
            % 3.1 compare old and new models
            f1 = figure('units','normalized','outerposition',[0 0 1 1]);
            for j = 1:length(rxns);
                subplot(2,5,j);
                bar(d, cell2mat(model1.(rxns{j})), w1, 'FaceColor',[0.2 0.2 1]);
                hold on
                bar(d, cell2mat(model2.(rxns{j})), w2, 'FaceColor',[1 0 0]);
                title(rxns{j});
                set(gca,'FontSize',16);
                xlim([-50,710])
                set(gca,'XTick',[50 250 450 650]);
                hold off
            end
            legend(sets{1}, sets{2},...
                'Orientation', 'horizontal', ...
                'Position',[0.4 0.50 0.22 0.041], 'FontSize',20, 'box', 'off');
            xlabel({'ATPM (\mumol/gDW/h)'},...
                'FontSize',30,'Position', [-1608 -0.539 0]);
            ylabel({'Contribution to total ATP production rate (%)'},...
                'FontSize',30);%, 'Position', [-4167.26 6.863 0]
            %%%%%%%%%
            %%%%%%%%%
            %%%%%%%%%
            %%% 3.2 stacked bar chart between old and new models
            w1 = 0.4;
            str = {'%'};
            f2= figure('units','normalized','outerposition',[0 0 1 1]);
            
           ydata = [cell2mat(model1.ATPS4mi), cell2mat(model1.glycolysis), (cell2mat(model1.r0408) ...
                + cell2mat(model1.RE0124C)), cell2mat(model1.SUCOASm),...
                (cell2mat(model1.NDPK6) + cell2mat(model1.NDPK2) + cell2mat(model1.UMPK) + cell2mat(model1.URIDK3)),...
                cell2mat(model1.NMNATr)];
            for j = 1:size(ydata,1);
                odata(j,1) = 100 - sum(ydata(j,:));
            end
            ydata(:,7) = odata;
            % Define a custom color palette with more colors
            customColors = [0 0.4470 0.7410;   % Blue
                0.8500 0.3250 0.0980;   % Orange
                0.9290 0.6940 0.1250;   % Yellow
                0.4940 0.1840 0.5560;   % Purple
                0.4660 0.6740 0.1880;   % Green
                0.3010 0.7450 0.9330;   % Cyan
                % Brown 0.8392 0.6941 0.5569;   
                0.6350 0.0780 0.1840];   % Red
            
            bars=bar(d+5, ydata, w1, 'stacked');
            colormap(customColors);
            for x = 1:numel(bars)
                set(bars(x),'Facecolor',customColors(x,:));
            end
            for m=1:size(ydata,1);
                xpos = d(m)+5;
                for n=1:size(ydata,2);
                    if ydata(m,n)>3; %so that the labels don't overlap with other labels
                        labels_stacked=num2str(ydata(m,n),'%.1f');
                        label = strcat(labels_stacked, '%');
                        hText = text(xpos, sum(ydata(m,1:n),2), label);
                        set(hText, 'VerticalAlignment','top', 'HorizontalAlignment',...
                            'center','FontSize',8, 'Color','k','FontWeight', 'Bold');
                    end
                end
            end
            hold on
             ydataPD = [cell2mat(model2.ATPS4mi), ...
                cell2mat(model2.glycolysis),(cell2mat(model2.r0408) ...
                + cell2mat(model2.RE0124C)), cell2mat(model2.SUCOASm),...
                (cell2mat(model2.NDPK6) + cell2mat(model2.NDPK2) + cell2mat(model2.UMPK) + cell2mat(model2.URIDK3)),...
                cell2mat(model2.NMNATr)];
            for m = 1:size(ydataPD,1);
                odataPD(m,1) = 100 - sum(ydataPD(m,:));
            end
            ydataPD(:,7) = odataPD;
            
            bars=bar(d+25, ydataPD, w1, 'stacked');
            colormap(customColors);
            for x = 1:numel(bars)
                set(bars(x),'Facecolor',customColors(x,:));
            end
            for m=1:size(ydataPD,1);
                xpos = d(m)+25;
                for n=1:size(ydataPD,2);
                    if ydataPD(m,n)>3; %so that the labels don't overlap with other labels
                        labels_stacked=num2str(ydataPD(m,n),'%.1f');
                        label = strcat(labels_stacked, '%');
                        hText = text(xpos, sum(ydataPD(m,1:n),2), label);
                        set(hText, 'VerticalAlignment','top', 'HorizontalAlignment',...
                            'center','FontSize',8, 'Color','k', 'FontWeight', 'Bold');
                    end
                end
            end
            set(gca,'FontSize',16, 'XTick',...
                [16 66 116 166 216 266 316 366 416 466 516 566 616 666], ...
                'XTickLabel', {...
                '      old  New','      old  New','      old  New','      old  New',...
                '      old  New','      old  New','      old  New','      old  New',...
                '      old  New','      old  New','      old  New','      old  New',...
                '      old  New','      old  New'});
            lgd =legend({'oxidative phosphorylation (ATPS4mi)', 'glycolysis (PGK & PYK)', ...
                'pentose phosphate pathway (r0408 & RE0124C)', ...
                'citric acid cycle (SUCOASm)',...
                'nucleotide interconversion (NDPK6 & NDPK2 & UMPK & URIDK3)', 'NAD metabolism (NMNATr)',...
                'other bioenergetic reactions'}, 'FontSize',10, 'box', 'off');
            %             lgd.Location='best';
            lgd.NumColumns=4;
            title([sets{1} ' vs ' sets{2}],'FontSize',30);
            xlabel({'ATPM (\mumol/gDW/h)'},'FontSize',30, ...
                'Position', [334 -7.5 0])
            ylabel({'Contribution to total ATP production rate (%)'},'FontSize',30);
            xlim([-10,710]);
            ylim([0,108]);
            set(gca,'XTickLabelRotation',0);
            %set(gca,'XTick',[])
            figure1 = f2;
            annotate
            caxis([0 7])
        end
    end
end

