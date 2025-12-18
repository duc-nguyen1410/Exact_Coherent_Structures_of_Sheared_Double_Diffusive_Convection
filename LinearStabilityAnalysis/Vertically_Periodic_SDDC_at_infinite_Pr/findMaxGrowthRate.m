function maxGR = findMaxGrowthRate(Ri,Pe,Lambda,Pr,tau)
    % adaptive resolutions based on Pe and Ri
    optN = 64;
    % using reduced SDDC
    % if 1e3<=Pe && Ri<2
    %     optN = 64; 
    % elseif Pe<10 && Ri>10
    %     optN = 80;
    % end
    % using full SDDC
    if 1e2<=Pe && Pe<1e3 && Ri<2
        optN = 100;
    elseif 1e3<=Pe && Ri<2
        optN = 200; 
    end

    % adaptive wavenumber range
    kx_list=[linspace(0.001,50/Pe,50) linspace(2,4,50)]';
    if Pe<25
        kx_list=[linspace(0.001,2,50)]';
    end
    sizeOfkx = size(kx_list,1);
    ky=0;

    maxGR = -inf;
    for kx_index=1:sizeOfkx
        kx=kx_list(kx_index);
        % compute eigenvalues
        % [~,eig_val] = eig_reducedModel(Ri,Pe,Lambda,Pr,tau,kx,ky,optN); % use Radko2016's reduced equations
        % [~,eig_val] = eig_Radko2016(Ri,Pe,Lambda,Pr,tau,kx,ky,optN); % use Radko2016's full equations
        [~,eig_val] = eig_fullSDDC(Ri,Pe,Lambda,Pr,tau,kx,ky,optN); % use full SDDC equations using Ra
        % [~,eig_val] = eig_reducedSDDC(Ri,Pe,Lambda,Pr,tau,kx,ky,optN); % use reduced SDDC equations using Ra
        eig_val(find(real(eig_val)>10^5))=-Inf;
        % compute growth rate
        if max(real(diag(eig_val))) > maxGR
            maxGR = max(real(diag(eig_val)));
        end
    end
end