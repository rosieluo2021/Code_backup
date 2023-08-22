function  plotIntersectRxns(ATPcontribution)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

new = {'ATPS4minew'
    'CYOOm2inew'
    'CYOOm3inew'
    'CYOR_u10minew'
    'NADH2_u10minew'
    'EX_cys_L(e)'};

original = {'ATPS4mi'
    'CYOOm2i'
    'CYOOm3i'
    'CYOR_u10mi'
    'NADH2_u10mi'};

allmodels=ATPcontribution;

types=fieldnames(allmodels);
for i=1:length(types)
    models=fieldnames(allmodels.(types{i}));
    models=models(~contains(models,'constrain'));% ignore unconstrained models
    if length(models)==2
        for m=1:length(new)
            A = find(ismember(allmodels.(types{i}).(models{1}).metRs(:,1),new(m,1)));
            if ~isempty(A)
                allmodels.(types{i}).(models{1}).metRs(A,1)= original(m,1);
            end
            B = find(ismember(allmodels.(types{i}).(models{2}).metRs(:,1),new(m,1)));
            if ~isempty(B)
                allmodels.(types{i}).(models{2}).metRs(B,1)= original(m,1);
            end
            clear A B
        end
        
        interRxns=intersect(allmodels.(types{i}).(models{1}).metRs(:,1),allmodels.(types{i}).(models{2}).metRs(:,1));
        fluxValue=zeros(length(interRxns),2);
        for j=1:length(interRxns)
            bool1=ismember(allmodels.(types{i}).(models{1}).metRs(:,1),interRxns{j});
            bool2=ismember(allmodels.(types{i}).(models{2}).metRs(:,1),interRxns{j});
            fluxValue(j,1)=allmodels.(types{i}).(models{1}).metRs{(bool1),2};
            fluxValue(j,2)=allmodels.(types{i}).(models{2}).metRs{(bool2),2};
        end
        % plot flux of the intersect rxns for two models
        c = string(interRxns(:,1));
            m= size(c,1);
            figure
            b=bar(fluxValue);%,'FaceColor' , [0 0.4470 0.7410]
            set(gca,'XTick',[1:m],'xticklabel',c)
            set(gca,'XTickLabelRotation',20, 'FontSize',12)
            ylabel('Flux(umol/gDW/h)', 'FontSize',14,'FontWeight','bold')
            l0=char(string(length(interRxns)));
            title(['overlapped ' l0 ' ATP contributing rxns'],'FontSize',12);
            l1=char(string(length(allmodels.(types{i}).(models{1}).metRs(:,1))));
            l2=char(string(length(allmodels.(types{i}).(models{2}).metRs(:,1))));
            legend([models{1} '(total ' l1 ' ATP rxns' ],[models{2} '(total ' l2 ' ATP rxns' ]);
            h = gca;
            h.YAxis.FontWeight = 'bold';
            h.XAxis.FontWeight = 'bold';
   
            for i=1:length(b)
                text(b(i).XEndPoints,b(i).YEndPoints+0.5,string(fluxValue(:,i)), ...
                          'VerticalAlignment','bottom','horizontalalign','center')
            end
            legend('boxoff')
    end
end
end

