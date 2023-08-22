function plotyy_eFBA(model,objective, max_C, N)
% This function plot variables of entropicFBA in a plotyy
%  Flux through biomass reaction (solution.v)) as a function of C_value on axis1 and
%  and non- linear part of objective function as a function od C_value on axis2
%
% inputs: 
%           model: a metabolic model that contain required fields to
%                  perform entropicFluxBalanceAnalysis
%           max_C : An estimation for the maximum value of C_value ( for
%                   larger number the variables do not change)
%           N : The number of C_value (since the model my be infeasible for some value of C 
%               the plotted figure will have fewer numbers in x axis

%  BY Samira Ranjbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FBAsolution = optimizeCbModel(model,'max');

param.solver = 'mosek';
param.printLevel = 0;
if strcmp(param.solver,'mosek')
    %set default mosek parameters for this type of problem
    param = mosekParamSetEFBA(param);
end

rng(0,'twister');
m = randi([1 max_C],1,N);
j = 1;

for i = 1:length(m)
    if any(model.c ~= 0)
    model.c(model.c ~= 0) = m(i);%?
    [solution,~] = entropicFluxBalanceAnalysis(model,param);
    if solution.stat == 1
        C(j)= m(i);
        j=j+1;
    end
    else
        error('no objective function.')
    end
end

C = sort(C);
solution_vals = ones(length (C),4);
param.printLevel = 0;

for nr = 1:length(C)
    model.c(model.c ~= 0) = C(nr);
    [solution,~] = entropicFluxBalanceAnalysis(model,param);
    solution_vals(nr,1) = solution.v(ismember(model.rxns, objective));
    solution_vals(nr,2) = solution.obj;
    solution_vals(nr,3) = solution.objEntropy;
    solution_vals(nr,4) = solution.objLinear; 
end


%%%%%%%%%%      solution.v  & solution.obj      %%%%%%%%%%%%%%%%%%%
[ax, h1, h2] = plotyy(C, solution_vals(:,1), C, solution_vals(:,2));

ylabel(ax(1), 'Objective flux (umol/gDw)', "FontSize",12, "FontWeight","bold");
ylabel(ax(2), 'Total Objective function', "FontSize",12, "FontWeight","bold");
xlabel('C-value', "FontSize",12, "FontWeight","bold")

grid on
grid minor

xlim([0 max_C])
set(h1, 'LineStyle', '-', 'Marker', 'o', 'Color', 'b', 'LineWidth',2);
set(h2, 'LineStyle', '--', 'Marker', 's', 'Color', 'r', 'LineWidth',2);

title('Effect of increasing c on predicted biomass')
set(gca, 'FontSize',12)

%hold on
%h3 = yline(FBAsolution.f,'--','MaxBiomasss','LineWidth',3);
%h3.LabelHorizontalAlignment = 'left';
%h3.Color=[0 1 0];
%hold off
hold on
h3 = FBAsolution.f*ones(size(1:max_C));
plot(1:max_C, h3,'Color', 'green','LineWidth',2)
hold off

legend('Objective-flux(eFBA)', 'Max objective flux(FBA)',' Objective function' , 'FontWeight','bold')
legend('Location', 'best') 


%%%%%%%%%%      solution.v  & solution.objEntropy      %%%%%%%%%%%%%%%%%%%
figure
[ax, h1, h2] = plotyy(C, solution_vals(:,1), C, solution_vals(:,3));

ylabel(ax(1), 'Objective flux (umol/gDw)', "FontSize",12, "FontWeight","bold");
ylabel(ax(2), 'Entropic part of objective function', "FontSize",12, "FontWeight","bold");
xlabel('C-value', "FontSize",12, "FontWeight","bold")

grid on
grid minor

xlim([0 max_C])
set(h1, 'LineStyle', '-', 'Marker', 'o', 'Color', 'b', 'LineWidth',2);
set(h2, 'LineStyle', '--', 'Marker', 's', 'Color', 'r', 'LineWidth',2);

title('Effect of increasing c on predicted biomass')
set(gca, 'FontSize',12)

hold on
h3 = FBAsolution.f*ones(size(1:max_C));
plot(1:max_C, h3,'Color', 'green','LineWidth',2)
hold off

legend('Objective-flux(eFBA)', 'Max objective flux(FBA)','Entropic part of objective function' , 'FontWeight','bold')
legend('Location', 'best') 

%%%%%%%%%%      solution.v  & solution.objLinear      %%%%%%%%%%%%%%%%%%%
figure
[ax, h1, h2] = plotyy(C, solution_vals(:,1), C, solution_vals(:,4));

ylabel(ax(1), 'Objective flux (umol/gDw)', "FontSize",12, "FontWeight","bold");
ylabel(ax(2), 'Linear part of objective function', "FontSize",12, "FontWeight","bold");
xlabel('C-value', "FontSize",12, "FontWeight","bold")

grid on
grid minor

xlim([0 max_C])
set(h1, 'LineStyle', '-', 'Marker', 'o', 'Color', 'b', 'LineWidth',2);
set(h2, 'LineStyle', '--', 'Marker', 's', 'Color', 'r', 'LineWidth',2);

title('Effect of increasing c on predicted biomass')
set(gca, 'FontSize',12)

hold on
h3 = FBAsolution.f*ones(size(1:max_C));
plot(1:max_C, h3,'Color', 'green','LineWidth',2)
hold off

legend('Objective-flux(eFBA)', 'Max objective flux(FBA)',' Linear part of objective function' , 'FontWeight','bold')
legend('Location', 'best')
end

