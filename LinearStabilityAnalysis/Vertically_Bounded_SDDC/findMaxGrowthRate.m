function maxGR = findMaxGrowthRate(Ra,Ri,Pr,Lambda,tau)
    % adaptive resolutions 
    optN = 128;
    % adaptive wavenumber ranges 
    opt_max_kx = 1.0;
    kx_list=linspace(0.001,opt_max_kx,2)';
    sizeOfkx = size(kx_list,1);
    ky=0;

    maxGR = -inf;
    for kx_index=1:sizeOfkx
        kx=kx_list(kx_index);
        % compute eigenvalues
        [~,eig_val] = eig_SDDC(Ra,Ri,Pr,Lambda,tau,kx,ky,optN);
        eig_val(find(real(eig_val)>10^2))=-Inf;
        % compute growth rate
        if max(real(diag(eig_val))) > maxGR
            maxGR = max(real(diag(eig_val)));
        end
    end
end