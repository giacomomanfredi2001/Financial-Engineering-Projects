function G = simulation_G(F0, TTM, k, nSim, flag)
% Simulation of the G for the MC computation
% 
% INPUT:
% F0:               value of the Forward at initial sim time
% TTM:              time to maturity till the end of the sim
% k:                volvol of the NMVM model
% nSim:             number of simulations
% flag:             [1 NIG, 2: VG]

    if flag == 1
        G = random('InverseGaussian', 1, TTM/k, [nSim, length(F0)]);
    else
        G = random('Gamma', TTM/k, k/TTM, [nSim, length(F0)]);
    end
    
    % Verification of the moments of 1st and 2nd order
    meanG = mean(mean(G));
    varG = mean((std(G)).^2);

end % function simulation_G