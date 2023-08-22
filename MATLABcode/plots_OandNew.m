function plots_OandNew(demandATP,f)
%UNTITLED3 Summary of this function goes here
%   change demandATP structure
if ~exist('f','var')
    f='barchart';
end

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
    
    %%%%%%
    
    %The reactions and X-axis (ATP demand) are the same for all models
    rxns = {'ATPS4mi', 'glycolysis', 'PGK', 'PYK', 'r0408', ...
        'RE0124C', 'SUCOASm', 'NDPK6', 'NDPK2', 'NMNATr'};
    if any(contains(models,'PD'))
        % Old and new
            type.index1=~contains(models,'PD');% old ane new control models
            type.index2=contains(models,'PD');% old ane new PD models
        
        for k=1:length(fieldnames(type))
            name=fieldnames(type);
            sets=models(type.(name{k}));
            model1=demandATP.(types{i}).(sets{1});
            model2=demandATP.(types{i}).(sets{2});
            
            customColors1 = [0 0.4470 0.7410;   % Blue
                0.8500 0.3250 0.0980;   % Orange
                0.9290 0.6940 0.1250;   % Yellow
                0.4940 0.1840 0.5560;   % Purple
                0.4660 0.6740 0.1880;   % Green
                0.8392 0.6941 0.5569;    % Brown
                0.6350 0.0780 0.1840];   % Red
            % Cyan 0.3010 0.7450 0.9330;
            
            customColors2 = [0.3010 0.7450 0.9330;   % Cyan(light blue)
                1     0.41  0.16;   % light Orange
                0.92  0.79  0.48;   % light Yellow
                0.72  0.27  1;   % light Purple
                0.39  0.83  0.07;   % light Green
                0.81  0.58  0.37;    % light Brown
                0.93  0.52 0.6];   % light Red
            
            Pathways={'oxidative phosphorylation (ATPS4mi)', 'glycolysis (PGK & PYK)', ...
                'pentose phosphate pathway (r0408 & RE0124C)', ...
                'citric acid cycle (SUCOASm)',...
                'nucleotide interconversion (NDPK6 & NDPK2 & UMPK & URIDK3)', 'NAD metabolism (NMNATr)',...
                'other bioenergetic reactions'};
            
            % y data for each model
            if ~isempty(cell2mat(model1.UMPK))
            ydata = [cell2mat(model1.ATPS4mi), cell2mat(model1.glycolysis), (cell2mat(model1.r0408) ...
                + cell2mat(model1.RE0124C)), cell2mat(model1.SUCOASm),...
                (cell2mat(model1.NDPK6) + cell2mat(model1.NDPK2) + cell2mat(model1.UMPK) + cell2mat(model1.URIDK3)),...
                cell2mat(model1.NMNATr)];
            else
            ydata = [cell2mat(model1.ATPS4mi), cell2mat(model1.glycolysis), (cell2mat(model1.r0408) ...
                + cell2mat(model1.RE0124C)), cell2mat(model1.SUCOASm),...
                (cell2mat(model1.NDPK6) + cell2mat(model1.NDPK2) + cell2mat(model1.URIDK3)),...
                cell2mat(model1.NMNATr)];
            end
            for j = 1:size(ydata,1);
                odata(j,1) = 100 - sum(ydata(j,:));
            end
            ydata(:,7) = odata;
            
            if ~isempty(cell2mat(model2.UMPK))
            ydataPD = [cell2mat(model2.ATPS4mi), ...
                cell2mat(model2.glycolysis),(cell2mat(model2.r0408) ...
                + cell2mat(model2.RE0124C)), cell2mat(model2.SUCOASm),...
                (cell2mat(model2.NDPK6) + cell2mat(model2.NDPK2) + cell2mat(model2.UMPK) + cell2mat(model2.URIDK3)),...
                cell2mat(model2.NMNATr)];
            else
                ydataPD = [cell2mat(model2.ATPS4mi), ...
                cell2mat(model2.glycolysis),(cell2mat(model2.r0408) ...
                + cell2mat(model2.RE0124C)), cell2mat(model2.SUCOASm),...
                (cell2mat(model2.NDPK6) + cell2mat(model2.NDPK2) + cell2mat(model2.URIDK3)),...
                cell2mat(model2.NMNATr)];
            end
            for m = 1:size(ydataPD,1);
                odataPD(m,1) = 100 - sum(ydataPD(m,:));
            end
            ydataPD(:,7) = odataPD;
            
            %%%%%%
            %%%%%%
            %%%%%% 1 split the stack bars to independent bars
            if f=='barchart'
                w1 = 0.4;
                str = {'%'};
                col=size(ydata,2);
                f2= figure('units','normalized','outerposition',[0 0 1 1]);
                for j = 1:col
                    subplot(col,1,j);
                    % add model1 plot
                    ydata_new=ydata(:,j);
                    bars=bar(d+5, ydata_new,w1, 'stacked','FaceColor', 'flat');
                    colormap(customColors1);
                    set(bars,'Facecolor',customColors1(j,:));
                    set(gca, 'box','off')
                    for m=1:size(ydata_new,1);
                        xpos = d(m)+5;
                        label=num2str(ydata_new(m),'%.1f');
                        label = strcat(label, '%');
                        hText = text(xpos, ydata_new(m), label);
                        set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                            'center','FontSize',8, 'Color','k', 'FontWeight', 'Bold');
                    end
                    
                    hold on
                    % add model2 plot
                    ydataPD_new=ydataPD(:,j);
                    bars=bar(d+25, ydataPD_new,w1, 'stacked','FaceColor', 'flat');
                    colormap(customColors2);
                    set(bars,'Facecolor',customColors2(j,:));
                    ylim([0,max(ydataPD_new)+10]);
                    set(gca, 'box','off')
                    for m=1:size(ydataPD_new,1);
                        xpos = d(m)+25;
                        label=num2str(ydataPD_new(m),'%.1f');
                        label = strcat(label, '%');
                        hText = text(xpos, ydataPD_new(m), label);
                        set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                            'center','FontSize',8, 'Color','k', 'FontWeight', 'Bold');
                    end
                    ylim([0,max([ydata_new; ydataPD_new])+5]);
                    xlim([-10,710]);
                    set(gca,'XTickLabelRotation',0);
                    title(Pathways{j},'FontSize',10);
                    set(gca,'FontSize',12, 'XTick',...
                        [16 66 116 166 216 266 316 366 416 466 516 566 616 666], ...
                        'XTickLabel', {...
                        '        old     new','        old     new','        old     new','        old     new',...
                        '        old     new','        old     new','        old     new','        old     new',...
                        '        old     new','        old     new','        old     new','        old     new',...
                        '        old     new','        old     new'});
                    %                     lgd=legend({sets{1},sets{2}}, 'FontSize',10, 'box', 'off');
                    %                     lgd.NumColumns=2;
                    if j~=7
                        set(gca,'XTick',[]);
                    end
                    if j==7
                        xlabel({'ATPM (umol/gDW/h)'},'Position', [345 -8.5 0],'FontSize',20)
                        figure1 = f2;
                        annotate
                    end
                    % add annotation
                    %                 annotation('textbox',...
                    %                     [0.45 0.90 0.2 0.08],...
                    %                     'String',{[sets{1} ' vs ' sets{2}]},...
                    %                     'FontSize',20,'FitBoxToText','on','LineStyle','none');
                end
                    % Give common xlabel, ylabel and title to your figure
                    han=axes(f2,'visible','off');
                    han.Title.Visible='on';
                    han.XLabel.Visible='on';
                    han.YLabel.Visible='on';
                    ylabel(han,'Contribution to total ATP production rate (%)','FontSize',20);
                    title(han,[sets{1} ' vs ' sets{2}],'FontSize',20);
            end
            %%%%%%%%
            %%%%%%%%
            %%%%%%%% 2. pie chart for each model
            if f=='piechart'
                col=size(ydata,1);
                f2= figure('units','normalized','outerposition',[0 0 1 1]);
                n=1;
                %explode = [1 1 1 1 1 1 1];
                for j = 1:col
                    subplot(col,2,n);
                    ydata_new=ydata(j,:);
                    pie(ydata_new);
                    n=n+1;
                    subplot(col,2,n);
                    ydataPD_new=ydataPD(j,:);
                    pie(ydataPD_new);
                    %title({['ATPM=' char(string(d(j))) '(umol/gDW/h)']},'FontSize',10);
                    n=n+1;
                end
            end
            
        end
    end
    
    if ~any(contains(models,'PD')) & any(contains(models,'ASTRO'))
        % Old and new
            sets=fieldnames(demandATP.(types{i}));
            model1=demandATP.(types{i}).(sets{1});
            model2=demandATP.(types{i}).(sets{2});
            
            customColors1 = [0 0.4470 0.7410;   % Blue
                0.8500 0.3250 0.0980;   % Orange
                0.9290 0.6940 0.1250;   % Yellow
                0.4940 0.1840 0.5560;   % Purple
                0.4660 0.6740 0.1880;   % Green
                0.8392 0.6941 0.5569;    % Brown
                0.6350 0.0780 0.1840];   % Red
            % Cyan 0.3010 0.7450 0.9330;
            
            customColors2 = [0.3010 0.7450 0.9330;   % Cyan(light blue)
                1     0.41  0.16;   % light Orange
                0.92  0.79  0.48;   % light Yellow
                0.72  0.27  1;   % light Purple
                0.39  0.83  0.07;   % light Green
                0.81  0.58  0.37;    % light Brown
                0.93  0.52 0.6];   % light Red
            
            Pathways={'oxidative phosphorylation (ATPS4mi)', 'glycolysis (PGK & PYK)', ...
                'pentose phosphate pathway (r0408 & RE0124C)', ...
                'citric acid cycle (SUCOASm)',...
                'nucleotide interconversion (NDPK6 & NDPK2 & UMPK & URIDK3)', 'NAD metabolism (NMNATr)',...
                'other bioenergetic reactions'};
            
            % y data for each model
            if ~isempty(cell2mat(model1.UMPK))
            ydata = [cell2mat(model1.ATPS4mi), cell2mat(model1.glycolysis), (cell2mat(model1.r0408) ...
                + cell2mat(model1.RE0124C)), cell2mat(model1.SUCOASm),...
                (cell2mat(model1.NDPK6) + cell2mat(model1.NDPK2) + cell2mat(model1.UMPK) + cell2mat(model1.URIDK3)),...
                cell2mat(model1.NMNATr)];
            else
            ydata = [cell2mat(model1.ATPS4mi), cell2mat(model1.glycolysis), (cell2mat(model1.r0408) ...
                + cell2mat(model1.RE0124C)), cell2mat(model1.SUCOASm),...
                (cell2mat(model1.NDPK6) + cell2mat(model1.NDPK2) + cell2mat(model1.URIDK3)),...
                cell2mat(model1.NMNATr)];
            end
            for j = 1:size(ydata,1);
                odata(j,1) = 100 - sum(ydata(j,:));
            end
            ydata(:,7) = odata;
            
            if ~isempty(cell2mat(model2.UMPK))
            ydataPD = [cell2mat(model2.ATPS4mi), ...
                cell2mat(model2.glycolysis),(cell2mat(model2.r0408) ...
                + cell2mat(model2.RE0124C)), cell2mat(model2.SUCOASm),...
                (cell2mat(model2.NDPK6) + cell2mat(model2.NDPK2) + cell2mat(model2.UMPK) + cell2mat(model2.URIDK3)),...
                cell2mat(model2.NMNATr)];
            else
                ydataPD = [cell2mat(model2.ATPS4mi), ...
                cell2mat(model2.glycolysis),(cell2mat(model2.r0408) ...
                + cell2mat(model2.RE0124C)), cell2mat(model2.SUCOASm),...
                (cell2mat(model2.NDPK6) + cell2mat(model2.NDPK2) + cell2mat(model2.URIDK3)),...
                cell2mat(model2.NMNATr)];
            end
            for m = 1:size(ydataPD,1);
                odataPD(m,1) = 100 - sum(ydataPD(m,:));
            end
            ydataPD(:,7) = odataPD;
            
            %%%%%%
            %%%%%%
            %%%%%% 1 split the stack bars to independent bars
            if f=='barchart'
                w1 = 0.4;
                str = {'%'};
                col=size(ydata,2);
                f2= figure('units','normalized','outerposition',[0 0 1 1]);
                for j = 1:col
                    subplot(col,1,j);
                    % add model1 plot
                    ydata_new=ydata(:,j);
                    bars=bar(d+5, ydata_new,w1, 'stacked','FaceColor', 'flat');
                    colormap(customColors1);
                    set(bars,'Facecolor',customColors1(j,:));
                    set(gca, 'box','off')
                    for m=1:size(ydata_new,1);
                        xpos = d(m)+5;
                        label=num2str(ydata_new(m),'%.1f');
                        label = strcat(label, '%');
                        hText = text(xpos, ydata_new(m), label);
                        set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                            'center','FontSize',8, 'Color','k', 'FontWeight', 'Bold');
                    end
                    
                    hold on
                    % add model2 plot
                    ydataPD_new=ydataPD(:,j);
                    bars=bar(d+9, ydataPD_new,w1, 'stacked','FaceColor', 'flat');
                    colormap(customColors2);
                    set(bars,'Facecolor',customColors2(j,:));
                    ylim([0,max(ydataPD_new)+10]);
                    set(gca, 'box','off')
                    for m=1:size(ydataPD_new,1);
                        xpos = d(m)+9;
                        label=num2str(ydataPD_new(m),'%.1f');
                        label = strcat(label, '%');
                        hText = text(xpos, ydataPD_new(m), label);
                        set(hText, 'VerticalAlignment','bottom', 'HorizontalAlignment',...
                            'center','FontSize',8, 'Color','k', 'FontWeight', 'Bold');
                    end
                    ylim([0,max([ydata_new; ydataPD_new])+5]);
                    xlim([0,55]);
                    set(gca,'XTickLabelRotation',0);
                    title(Pathways{j},'FontSize',10);
                    set(gca,'FontSize',12, 'XTick',...
                        [7 17 27 37 47], ... 
                        'XTickLabel', {...
                        '             old            new',...
                        '             old            new',...
                        '             old            new',...
                        '             old            new',...
                        '             old            new'});

                    if j~=7
                        set(gca,'XTick',[]);
                    end
                    if j==7
                        xlabel({'ATPM (umol/gDW/h)'},'Position', [345 -8.5 0],'FontSize',20)
                        figure1 = f2;
                        annotate1
                    end
                  
                end
                    % Give common xlabel, ylabel and title to your figure
                    han=axes(f2,'visible','off');
                    han.Title.Visible='on';
                    han.XLabel.Visible='on';
                    han.YLabel.Visible='on';
                    ylabel(han,'Contribution to total ATP production rate (%)','FontSize',20);
                    title(han,[sets{1} ' vs ' sets{2}],'FontSize',20);
            end
    end
end





