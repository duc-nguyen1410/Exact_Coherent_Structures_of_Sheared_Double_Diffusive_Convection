clear all;
close all;
clc;
setfigure;

% define paramesters 
Ra= 1e7;
Ri = 1000; % Richardson number
Pr = 7; % Prandtl number
Lambda = 2.; % density ratio
Le = 100; % Lewis number = 1/tau
filename = ['growth_rate_kxky'];
kx_list= linspace(0.0001,1,10)'; % horizontal wavenumber in x direction
ky_list= linspace(0.0001,1,10)'; % horizontal wavenumber in y direction';
N=128;

GR = getGrowthRate_SDDC_kxky(Ra,Ri,Pr,Lambda,Le,kx_list,ky_list,N);
% size(GR)
writematrix(GR,[filename '.csv']);
[kx,ky,interpGR] = interp(kx_list,ky_list,GR,100,100);
% [maxmax,idx] = max(interpGR,[],"all");
% [Ramax,Rimax] = ind2sub(size(interpGR),idx);
% disp(['max value = ' num2str(maxmax) ' of growth rate at (Ra,Ri)=(' num2str(abs(Ra(Ramax))) ',' num2str(abs(Ri(Rimax))) ')']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <---- export data in vtk-format
% vtkwrite([filename '.vtk'],'structured_points','growth_rate', transpose(interpGR));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpGR(find(interpGR<1e-5))=NaN; % NaN values will be displayed on figure -> white color
% plot growth rate
f = figure;
% plot(kx_list,GR);
% xlabel('{\it{k_x}}');
pcolor(kx,ky,transpose(interpGR));shading interp;colormap();colorbar;
% set(gca, 'ZScale', 'log')
pbaspect([1 1 1])
% caxis([-0.1, 0.1])
xlabel('{\it{kx}}')
ylabel('{\it{ky}}',"Rotation",0)
savefigure(gca,'growth_rate_kxky_small.png');