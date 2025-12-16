clear all;
close all;
clc;
setfigure;

% define paramesters
Ra = 1e7;
Ri = 10.;
Rrho = 2.;
Pr = 7.;
Le = 100; 
Lambda = 2.0; % density ratio
kx = 1;
ky = 0;

[eig_vec,eig_val] = eig_SDDC(Ra,Ri,Pr,Lambda,Le,kx,ky,128);

eig_val=diag(eig_val);
eig_val(find(abs(real(eig_val))>10^1))=NaN;
fprintf('max of Re(eig_val) is %s\n', mat2str(max(real(eig_val))));
% writematrix(eig_vec,'eig_vec_tab.txt','Delimiter','tab')
% writematrix(eig_val,'eig_val_tab.txt','Delimiter','tab')

% plot eigenvalues
f = figure;
plot(real(eig_val),imag(eig_val),'o');
line([0 0], ylim,'Color','black');  %x-axis
line(xlim, [0 0],'Color','black');  %y-axis
xlabel('real') 
ylabel('imag',"Rotation",0) 
% title(['eigenvalues Ri=' num2str(Ri),' Pe=' num2str(Pe)])
% savefigure(gca,['eigenvalues_Ri=' num2str(Ri) '_Pe=' num2str(Pe) '.png']);