function GR = getGrowthRate_SDDC_RaRi(Ra_list,Ri_list,Pr,Lambda,tau)
    % This function uses Parallel Computing Toolbox to rise the computing
    % speed. Normally, a single core will be used when Parallel Computing 
    % Toolbox is not installed on the local machine.
    sizeOfRa = size(Ra_list,1);
    sizeOfRi = size(Ri_list,1);
    % create an array to store values of growth rate
    GR = zeros(sizeOfRa,sizeOfRi);
    if canUseParallelPool
        % Parallel Computing Toolbox is installed
        parfor Ra_index=1:sizeOfRa
            Ra=Ra_list(Ra_index);
            for Ri_index=1:sizeOfRi
                Ri=Ri_list(Ri_index);
                GR(Ra_index,Ri_index)=findMaxGrowthRate(Ra,Ri,Pr,Lambda,tau);
            end
        end
    else
        % Parallel Computing Toolbox is not installed
        for Ra_index=1:sizeOfRa
            Ra=Ra_list(Ra_index);
            for Ri_index=1:sizeOfRi
                Ri=Ri_list(Ri_index);
                GR(Ra_index,Ri_index)=findMaxGrowthRate(Ra,Ri,Pr,Lambda,tau);
            end
        end
    end
end