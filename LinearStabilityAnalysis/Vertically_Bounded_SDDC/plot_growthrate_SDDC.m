clc;
Ra_list = linspace(1,1e7,40)';
Ri_list = linspace(0.001,10,40)'; % Richardson number
GR = readmatrix("Linear_stability_analysis/growth_rate_Pr=7_Lambda=2.csv");

[Ra,Ri,interpGR] = interp(Ra_list,Ri_list,GR,100,100);
[maxmax,idx] = max(interpGR,[],"all");
[Ramax,Rimax] = ind2sub(size(interpGR),idx);
disp(['max value = ' num2str(maxmax) ' of growth rate at (Ra,Ri)=(' num2str(abs(Ra(Ramax))) ',' num2str(abs(Ri(Rimax))) ')']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <---- export data in vtk-format
vtkwrite([filename '.vtk'],'structured_points','growth_rate', transpose(interpGR));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot growth rate
f = figure;
pcolor(Ra,Ri,transpose(interpGR));shading interp;colormap(bwr);colorbar;
% set(gca, 'ZScale', 'log')
% set(gca, 'ColorScale', 'log')
% pbaspect([1 1 1])
% caxis([-0.25, 0.25])
xlabel('{\it{Ra}}')
ylabel('{\it{Ri}}',"Rotation",0)
savefigure(gca,[filename '.png']);