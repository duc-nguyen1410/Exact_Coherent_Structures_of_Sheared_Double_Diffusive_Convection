function GR = getGrowthRate_SDDC_kxky(Ra,Ri,Pr,Lambda,Le,kx_list,ky_list,N)
    % This function uses Parallel Computing Toolbox to rise the computing
    % speed. Normally, a single core will be used when Parallel Computing 
    % Toolbox is not installed on the local machine.
    sizeOfkx = size(kx_list,1);
    sizeOfky = size(ky_list,1);
    % create an array to store values of growth rate
    GR = zeros(sizeOfkx,1);
    if canUseParallelPool
        % Parallel Computing Toolbox is installed
        parfor kx_index=1:sizeOfkx
            kx = kx_list(kx_index);
            for ky_index=1:sizeOfky
                ky = ky_list(ky_index);
                % compute eigenvalues
                [~,eig_val] = eig_SDDC(Ra,Ri,Pr,Lambda,Le,kx,ky,N);
                eig_val(find(real(eig_val)>10^1))=-Inf;
                % compute growth rate
                GR(kx_index,ky_index)=max(real(diag(eig_val)));
            end
        end
    else
        % Parallel Computing Toolbox is not installed
        for kx_index=1:sizeOfkx
            kx=kx_list(kx_index);
            for ky_index=1:sizeOfky
                ky = ky_list(ky_index);
                % compute eigenvalues
                [~,eig_val] = eig_SDDC(Ra,Ri,Pr,Lambda,Le,kx,ky,N);
                eig_val(find(real(eig_val)>10^1))=-Inf;
                % compute growth rate
                GR(kx_index,ky_index)=max(real(diag(eig_val)));
            end
        end
    end
end