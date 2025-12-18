clear all;
close all;
clc;
setfigure;

% define paramesters 
Ra_list = linspace(1,1e7,100)';
Ri_list = linspace(0.001,10,100)'; % Richardson number
Pr = 7; % Prandtl number
Lambda = 2.; % density ratio
Le = 100; % Lewis number = 1/tau
filename = ['./2D/Linear_stability_analysis/growth_rate_Pr=' num2str(Pr),'_Lambda=' num2str(Lambda)];
kx=1;
ky=0;
N=128;

GR = getGrowthRate_SDDC_RaRi(Ra_list,Ri_list,Pr,Lambda,Le,kx,ky,N);
% size(GR)
writematrix(GR,[filename '.csv']);
[Ra,Ri,interpGR] = interp(Ra_list,Ri_list,GR,100,100);
[maxmax,idx] = max(interpGR,[],"all");
[Ramax,Rimax] = ind2sub(size(interpGR),idx);
disp(['max value = ' num2str(maxmax) ' of growth rate at (Ra,Ri)=(' num2str(abs(Ra(Ramax))) ',' num2str(abs(Ri(Rimax))) ')']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <---- export data in vtk-format
% vtkwrite([filename '.vtk'],'structured_points','growth_rate', transpose(interpGR));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpGR(find(interpGR<1e-5))=NaN; % NaN values will be displayed on figure -> white color
% plot growth rate
f = figure;
pcolor(Ra,Ri,transpose(interpGR));shading interp;colormap();colorbar;
% set(gca, 'ZScale', 'log')
% pbaspect([1 1 1])
% caxis([-0.1, 0.1])
xlabel('{\it{Ra}}')
ylabel('{\it{Ri}}',"Rotation",0)
savefigure(gca,[filename '.png']);