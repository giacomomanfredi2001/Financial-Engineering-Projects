function [price_call, Ft] = E3_montecarlo(strike, sigma, TTM, k, eta, alpha, discount, ...
    F0, nSim, flag, G)
% Computation of the Call Prices through the MC method
% 
% INPUT:
% strike:                     strike of stop
% sigma:                      volatility of the contract
% TTM:                        time to maturity
% k:                          volvol of the contract
% eta:                        simmetry of the contract
% alpha:                      parameter to define a particular model
% discount:                   vec discount at trigger_dates
% F0:                         initial forward value
% flag:                       [1: NIG, 2: VG]
% G:                          simulation of the G
% 
% OUTPUT:
% price_call:                 price of the call at maturity
% Ft:                         forward of the call at maturity
% 
% USES:
% function LaplaceTransform (TTM, k, alpha, eta, sigma)
 
    %% Creation of the moneyness grid
    
    x_grid = log(F0/strike);
    
    %% Simulation of the g

    g = randn(nSim, length(x_grid));

    %% Creation of ft - forward exponent dynamic

    mu = -log(LaplaceTransform(TTM, k, alpha, eta, sigma, flag));
    
    ft = mu - sigma.^2 .* (0.5 + eta).*TTM.*G + sigma.*sqrt(TTM.*G) .*g;

    %% Computation of the forward

    Ft = F0 * exp(ft);

    %% Computation of the call
    
    indicator = Ft >= ones(nSim, length(x_grid)).*strike;
    price_call = discount .* (Ft .* indicator - strike.*indicator);

end % function E3_montecarlo