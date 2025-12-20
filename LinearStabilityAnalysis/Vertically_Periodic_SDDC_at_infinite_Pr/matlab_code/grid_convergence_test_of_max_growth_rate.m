clear all;
close all;
clc;
setfigure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change this
Ri = 0.25; % Richardson number
Pe = 1; % Peclet number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change resolution here
N_list=[40,64,80,100]';
legend_list = "N="+num2str(N_list);
filename = ['./LinearStabilityAnalysis/Vertically_Periodic_SDDC_at_infinite_Pr/convergence_test_of_max_growth_rate_Ri=' num2str(Ri),'_Pe=' num2str(Pe)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters
Lambda = 2.; % density ratio
tau = 0.01; % diffusivity ratio
Pr = 10.;
% k_list=[linspace(0.001,50/Pe,50) linspace(2,10,50)]';
k_list=[linspace(0.001,2,50)]';
l_list=0;
f = figure;
for N_index=1:length(N_list)
    N = N_list(N_index);
    fprintf('N=%d\n',N);
    GR = growthrate(Ri,Pe,Lambda,Pr,tau,k_list,l_list,N);
    [maxmax,kmax] = max(GR,[],"all");
    plot(k_list,GR','Color', [0 0 0 N_index/length(N_list)],'LineWidth',2);hold on;
    % plot(k_list(kmax),maxmax,'o','Color', [0 0 0 N_index/length(N_list)]);hold on;
end
legend(legend_list);
% pbaspect([1 0.5 1])
xlabel('{\it{k}}')
ylabel('Growth rate')
savefigure(gca,[filename '.png']);