set(0, 'Defaultfigurecolor',[1 1 1])
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 20)
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 20)

n = 256;  % number of colors
r = [linspace(0,1,n/2), ones(1,n/2)]';
g = [linspace(0,1,n/2), linspace(1,0,n/2)]';
b = [ones(1,n/2), linspace(1,0,n/2)]';
bwr = [r g b];