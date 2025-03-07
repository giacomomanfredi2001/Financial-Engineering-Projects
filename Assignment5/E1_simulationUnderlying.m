function ES = E1_simulationUnderlying(monitoring_rates, dividend, initial_underlying, rho, sigma, interdelta,  nSim)
% Computation of the underlying using GBM
% 
%INPUT:
% monitoring_rates:                 monitoring rates of the contract
% dividend:                         vector of yield (ENEL, AXA)
% initial_underlying:               underlying at t=0 (ENEL, AXA)
% rho:                              correlation coefficients
% sigma:                            volatilities (ENEL, AXA)
% interdelta:                       delta in ACT365
% nSim:                             number of simulations

    %% Initialization of the parameters
    

    interdelta = ones(nSim, length(monitoring_rates), length(sigma)) .* interdelta';
    sigmav = ones(nSim, length(monitoring_rates), length(sigma));
    dividendsv = ones(nSim, length(monitoring_rates), length(dividend));
    monitoring_rates = ones(nSim, length(monitoring_rates), length(sigma)) .* monitoring_rates';

    for i=1:length(sigma)
        sigmav(:, :, i) = sigma(i);
    end

    sigma = sigmav;

    for i=1:length(dividend)
        dividendsv(:, :, i) = dividend(i);
    end

    dividend = dividendsv;

    %% Initialization of the simulation

    % Computation of the iid multivariate normal samples
    meanVector = zeros(size(sigma, 3), 1);
    correlationMatrix = rho*ones(size(sigma, 3)) + (1-rho)*diag(ones(size(sigma, 3), 1));
    
    samples = zeros(nSim, size(monitoring_rates, 2), size(sigma, 3));

    for i = 1:size(monitoring_rates, 2)
        samples(:, i, :) = mvnrnd(meanVector, correlationMatrix, nSim);
    end
    
    % Initialization of the matrix
    ES = zeros(nSim, size(monitoring_rates, 2)+1, size(sigma, 3));

    % [ENEL0 ENEL1 ENEL2 ENEL3 ENEL4] ... [AXA0 AXA1 AXA2 AXA3 AXA4] on
    % different level (3rd dimension)

    %% Simulation od the underlying

    % Initialization of the ES
    ES(:, 1, :) = ones(size(ES, 1), size(ES, 3)) .* initial_underlying';

    % Fulfillment of the ES 
    for i = 2:size(monitoring_rates, 2)+1

        a = (monitoring_rates(:, i-1, :) - dividend(:, i-1, :) - 0.5.*sigma(:, i-1, :).^2).*interdelta(:, i-1, :);
        b = sigma(:, i-1, :).*sqrt(interdelta(:, i-1, :)).*samples(:, i-1, :);
        c = exp(a + b);
                                  
        ES(:, i, :) = ES(:, i-1, :) .* c;
    end

end % function E1_simulationUnderlying