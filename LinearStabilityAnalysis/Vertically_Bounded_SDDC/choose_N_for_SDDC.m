% define paramesters
Ra = 1e7;
Ri = 1000.0; % Richardson number
tau = 0.01; % Peclet number
Lambda = 2.0; % density ratio
Pr = 7.0;
kx=0.00001;
ky=0;
N = [16 32 64 96 128 256]';

for n_id = 1:length(N)
    n = N(n_id);
    [~,eig_val] = eig_SDDC(Ra,Ri,Pr,Lambda,tau,kx,ky,n);
    eig_val(find(real(eig_val)>1))=-Inf;
    eigv_arr = diag(eig_val)';
    real_matrix = real(eigv_arr);
    [maxvalue,ind] = max(real_matrix,[],"all");
        
    fprintf("N = %d --> max growth rate = %.15f + %.15fi\n",n,real(eigv_arr(ind)),imag(eigv_arr(ind)));
end

% Ra=1e7, Ri = 1, Le = 1e2, Pr=7, Lambda=2, and kx=1, ky=0;  for 2D domain
% N = 16 --> max growth rate = 0.007737519947247 + 0.000000000000000i
% N = 32 --> max growth rate = 0.002904366354960 + -0.000000000000000i
% N = 64 --> max growth rate = -0.004225980379278 + 0.000000000000001i
% N = 96 --> max growth rate = -0.012329297799492 + -0.478722104452513i
% N = 128 --> max growth rate = -0.012329297790247 + 0.478722104534697i -- Select --
% N = 256 --> max growth rate = -0.012329297790005 + -0.478722104534401i
% N = 384 --> max growth rate = -0.012329297787122 + 0.478722104535462i
% N = 512 --> max growth rate = -0.012329297788483 + -0.478722104511622i

% Ra=1e7, Ri = 1, Le = 1e2, Pr=7, Lambda=2, and kx=1, ky=1; for 3D domain
% N = 16 --> max growth rate = -0.004800752519424 + 0.456013046620884i
% N = 32 --> max growth rate = -0.005418356336894 + 0.240353832832037i
% N = 64 --> max growth rate = -0.009718046261617 + 0.192417951679473i
% N = 96 --> max growth rate = -0.012430444912835 + 0.478849471711953i
% N = 128 --> max growth rate = -0.012430444898973 + 0.478849471792873i -- Select --
% N = 256 --> max growth rate = -0.012430444898525 + 0.478849471792165i
% N = 384 --> max growth rate = -0.012430444895893 + 0.478849471792804i