function [SFBA_rxns,overlapped] = SparseRxns(model,param)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

modelnames=fieldnames(model);
% Finds the minimal set of reactions subject to a LP objective
for i=1:length(modelnames)
    sparseFBAModel = model.(modelnames{i});
    [vSparse,sparseRxnBool,essentialRxnBool]  = sparseFBA(sparseFBAModel);

    minRxns = sparseFBAModel.rxns(find(sparseRxnBool));
    %minLB = sparseFBAModel.lb(findRxnIDs(sparseFBAModel, minRxns));
    %minUB = sparseFBAModel.ub(findRxnIDs(sparseFBAModel, minRxns));
    SFBA_rxns.(modelnames{i})=minRxns;
end

statistic=zeros(length(modelnames));
%calculate overlapped rxns
for i=1:length(modelnames)
    for j=1:length(modelnames)
    overlapped.(modelnames{i}).(['with' modelnames{j}])= SFBA_rxns.(modelnames{i})(ismember(SFBA_rxns.(modelnames{i}),SFBA_rxns.(modelnames{j})));
    statistic(i,j)=length(overlapped.(modelnames{i}).(['with' modelnames{j}]));
    end
end
%statistic=cell2table(statistic);
%statistic.Properties.RowNames=modelnames;
%statistic.Properties.VariableNames=modelnames;
%use proportion to create heat map
if isfield(param,'plot')
    if param.plot==1
        figure('units','normalized','outerposition',[0 0 1 1])
        for i=1:size(statistic,1)
            [max_a(i),index(i)]=max(statistic(i,:));
        end
        xa=repmat(max_a',[1 size(statistic,1)]);
        pro=round(statistic./xa*100,2);
        ax = subplot(1,1,1);
        h = imagesc(ax, pro);
        ax.TickLength(1) = 0;
        % Create heatmap's colormap
        n=256;
        cmap = [linspace(.9,0,n)', linspace(.9447,.447,n)', linspace(.9741,.741,n)'];
        colormap(ax, cmap);
        colorbar(ax)
        hold on
        %add text label (proportion (accuracy))
        label1=reshape(pro',[],1);
        % accuracy_mets=accuracy.mets{:,2:end};
        % label2=round(reshape(accuracy_mets,[],1),2);
        % labels=append(string(label1),'%','(',string(label2),')');
        try
            %labels(find(label1==100))=append(string(max_a'),'(',string(label2(find(label1==100))),')');
            labels=append(string(label1));
            labels(find(label1==100))=append(string(max_a'));
        catch ME
            disp(ME)
        end
        [xTxt, yTxt] = ndgrid(1:size(statistic,1), 1:size(statistic,1));
        th = text(xTxt(:), yTxt(:), labels(:), ...
            'VerticalAlignment', 'middle','HorizontalAlignment','Center');
        set(ax,'XTick',1:size(statistic,1),'YTick',1:size(statistic,1))
        xticklabels(modelnames)
        yticklabels(modelnames)
        set(gca,'XTickLabelRotation',0);
        %
        title('Overlapped minimal rxns for each pair of models')
        annotation('textbox',...
            [0.75 0.92 0.2 0.08],...
            'String',{'heatmap:','Colorbar = overlapped proportion(%)','textlabel = overlapped proportion(%)','Diagonal number = minimal rxn size of each model'},...
            'FontSize',10,'FitBoxToText','on','LineStyle','none');
    end
end

end

