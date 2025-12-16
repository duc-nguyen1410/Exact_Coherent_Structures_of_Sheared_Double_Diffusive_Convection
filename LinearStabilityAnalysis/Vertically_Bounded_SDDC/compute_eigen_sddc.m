clear all;
close all;
clc;
setfigure;

% define paramesters
Ra=1e7;
Lambda=2;%Lambda
Le=100;%=1/tau
Pr=7;
Ri=1;
kx=1;%Lx=2pi
ky=0;%Ly=0
N=128;

[eigenvector,eigenvalue] = eig_SDDC(Ra,Ri,Pr,Lambda,Le,kx,ky,N);

eigenvalue=diag(eigenvalue);
eigenvalue(find(real(eigenvalue)>1e2))=NaN;
fprintf('max of Re(eig_val) is %s\n', mat2str(max(real(eigenvalue))));

% plot eigenvalues
f = figure;
plot(real(eigenvalue),imag(eigenvalue),'o');
xlim([-1,0.1])
ylim([-0.5,0.5])
line([0 0], ylim,'Color','black');  %x-axis
line(xlim, [0 0],'Color','black');  %y-axis
xlabel('real') 
ylabel('imag',"Rotation",0) 
title(['eigenvalues Ra=' num2str(Ra),' Pr=' num2str(Pr)])
% savefigure(gca,['eigenvalues_Ri=' num2str(Ri) '_Pe=' num2str(Pe) '.png']);