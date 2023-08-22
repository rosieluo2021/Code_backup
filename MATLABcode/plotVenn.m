function plotVenn(compareExRxnsTable)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
type=fieldnames(compareExRxnsTable);
for i=1:length(type)
    figure('units','normalized','outerposition',[0 0 1 1]);
    files=fieldnames(compareExRxnsTable.(type{i}));
    subplotLayout = [5,4];
    ax = gobjects(fliplr(subplotLayout));
    a=1;
    for j=1:length(files)
        eachfile=compareExRxnsTable.(type{i}).(files{j});
        data=eachfile(2:5,2:4);
        for m=1:4
            ax(a)=subplot(length(files),4,a);
            %             A=data{m,1};
            %             B=data{m,2};
            %             C=data{m,3};
            %             datasets = [A B];
            %             labels={A,B,C};
            %             [H , S]=venn(datasets,C,'FaceColor',{'c','w'},'FaceAlpha',{1,0.6},'EdgeColor','w');
            %             title(eachfile{m+1,1},'FontSize',8)
            %             axis off
            %Now label each zone
%             for n = 1:3
%                 if S.ZoneCentroid(n,1)<0
%                     S.ZoneCentroid(n,1)=S.ZoneCentroid(n,1)*1.5;
%                 end
%                 text(S.ZoneCentroid(n,1), S.ZoneCentroid(n,2), string(labels(n)));
%             end
            A=num2str(data{m,1});
            B=num2str(data{m,2});
            C=num2str(data{m,3});
            modelnames=splitString(eachfile{m+1,1});
            modelnames=modelnames([1,3],:);
            venn_new(2,sets=modelnames,labels={A,B,'c','d',C});
            a=a+1;
        end
    end
    titles = {'Totol reactions', 'Ex/Sink/DM reactions', 'Mitochondrial reactions','Total genes','Mitochondrial genes'};
    %
    ax = ax';
    % Reduce vertical space of axes just a bit to make room for "titles"
    axPos = cell2mat(get(ax, 'Position'));
    axPos(:,4) = axPos(:,4).*.86; % reduces vert space to 96% of height
    set(ax, {'Position'}, mat2cell(axPos, ones(numel(ax),1), 4))
    
    % Get upper position of each row of axes, normalized coordinates
    axPos = cell2mat(get(ax(:,1), 'Position'));
    axUpperPos = sum(axPos(:,[2,4]),2);  %upper pos.
    % Get center position for 1st row (assumes all rows have same center)
    axPos = cell2mat(get(ax(1,[1,end]),'Position'));
    axCenterPos = mean([axPos(1,1), sum(axPos(2,[1,3]))]);
    
    % Set annotation for each row of subplots
    titleHandles = gobjects(numel(titles),1);
    for k = 1:numel(titles)
        titleHandles = annotation('textbox','String',titles{k}, ...
            'Position', [axCenterPos, axUpperPos(k), 0, 0], ...
            'HorizontalAlignment', 'center','VerticalAlignment','bottom',...
            'LineStyle','none','FitBoxToText','on', ...
            'FontWeight',ax(1).Title.FontWeight, ... % matches title property
            'FontSize', 13, ...    % matches title property
            'FontName', ax(1).Title.FontName);    % matches title property
            %'Color', ax(1).Title.Color);             % matches title property
    end

end

end

