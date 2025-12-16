clear all;
close all;
clc;
setfigure;

% define paramesters 
Ra= 1e7;
Ri = 1; % Richardson number
Pr = 7; % Prandtl number
Lambda = 2.; % density ratio
Le = 100; % Lewis number = 1/tau
filename = ['growth_rate_kx'];
kx_list= linspace(0.01,100,80)'; % horizontal wavenumber in x direction
ky=0;
N=196;

GR = getGrowthRate_SDDC_kx(Ra,Ri,Pr,Lambda,Le,kx_list,ky,N);
% size(GR)
% writematrix(GR,[filename '.csv']);
% [Ra,Ri,interpGR] = interp(Ra_list,Ri_list,GR,100,100);
% [maxmax,idx] = max(interpGR,[],"all");
% [Ramax,Rimax] = ind2sub(size(interpGR),idx);
% disp(['max value = ' num2str(maxmax) ' of growth rate at (Ra,Ri)=(' num2str(abs(Ra(Ramax))) ',' num2str(abs(Ri(Rimax))) ')']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <---- export data in vtk-format
% vtkwrite([filename '.vtk'],'structured_points','growth_rate', transpose(interpGR));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpGR(find(interpGR<1e-5))=NaN; % NaN values will be displayed on figure -> white color
% plot growth rate
f = figure;
plot(kx_list,GR);
xlabel('{\it{k_x}}');
% pcolor(Ra,Ri,transpose(interpGR));shading interp;colormap();colorbar;
% set(gca, 'ZScale', 'log')
% pbaspect([1 1 1])
% caxis([-0.1, 0.1])
% xlabel('{\it{Ra}}')
% ylabel('{\it{Ri}}',"Rotation",0)
savefigure(gca,'growth_rate_kx.png');