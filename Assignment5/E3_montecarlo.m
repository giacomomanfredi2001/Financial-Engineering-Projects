function price_call = E3_montecarlo(x_grid, sigma, TTM, k, eta, alpha, discount, F0, nSim)
% Computation of the Call Prices through the MC method
% 
%INPUT:
% x_grid:                     logmoneyness grid of the quadrature
% sigma:                      volatility of the contract
% TTM:                        time to maturity
% k:                          volvol of the contract
% eta:                        simmetry of the contract
% alpha:                      parameter to define a particular model
% discount:                   discount at maturity
% F0:                         initial forward value
% 
%USES:
% function LaplaceTransform (TTM, k, alpha, eta, sigma)
 
    %% Creation of the strike grid
    
    strike_grid = (F0*exp(-x_grid))';
    
    %% Simulation of the G, g

    g = randn(nSim, length(x_grid));

    G = random('InverseGaussian', 1, TTM/k, [nSim, length(x_grid)]);
    
    % Verification of the moments of 1st and 2nd order
    meanG = mean(mean(G));
    varG = mean((std(G)).^2);

    %% Creation of ft - forward exponent dynamic

    mu = -log(LaplaceTransform(TTM, k, alpha, eta, sigma));
    
    ft = mu - sigma.^2 .* (0.5 + eta).*TTM.*G + sigma.*sqrt(TTM.*G) .*g;

    %% Computation of the forward

    Ft = F0 * exp(ft);

    %% Computation of the call
    
    indicator = Ft >= ones(nSim, length(x_grid)).*strike_grid;
    price_call = discount .* (Ft .* indicator - strike_grid.*indicator);

end % function E3_montecarlo