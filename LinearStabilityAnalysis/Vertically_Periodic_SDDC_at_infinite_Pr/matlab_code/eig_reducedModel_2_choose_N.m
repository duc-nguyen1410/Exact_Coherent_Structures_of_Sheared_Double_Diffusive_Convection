% define paramesters
Ri = 10; % Richardson number
Pe = 1e2; % Peclet number
Pr = 10.;
tau = 0.01; % diffusivity ratio
Lambda = 2.; % density ratio

% Ra = 1e5; % Rayleigh number
Ra = (4*pi*pi*Ri*Pe*Pe)/(Pr*(Lambda-1)); % Rayleigh number

kx=0.11;
ky=0;
N = [20 30 50 100 200]';

for n_id = 1:length(N)
    n = N(n_id);
    [~,eig_val] = eig_reducedModel_2(Ri,Ra,Pe,Lambda,Pr,tau,kx,ky,n);
    eig_val(find(real(eig_val)>1))=-Inf;
    eigv_arr = diag(eig_val)';
    real_matrix = real(eigv_arr);
    [maxvalue,ind] = max(real_matrix,[],"all");
        
    fprintf("N = %d --> max growth rate = %.15f + %.15fi\n",n,real(eigv_arr(ind)),imag(eigv_arr(ind)));
end

% Ri = 1.0, Pe = 1e2, and kx=0.27; 
% 