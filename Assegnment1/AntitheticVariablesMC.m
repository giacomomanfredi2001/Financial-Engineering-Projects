function [unbiasedStd, unbiasedStdAV] = AntitheticVariablesMC(F0, K, B, T, sigma, N)
%Computation of the Vega Greek for the hedging
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% N:     number of simulations
%
%OUTPUT
% unbiasedStd:        error without AV technique
% unbiasedStdAV:      error with AV technique

Nsim = 100;

% Computation of the values
optionPricesMC = [];
optionPricesMC_COMB = [];

for j = 1:Nsim
    g = randn(N, 1);      % Normal std distribution

    optionPriceVector = max(F0 .* exp(-0.5 * sigma^2 * T - sigma .* sqrt(T) .* g) - K, 0);
    optionPriceVectorAV = max(F0 .* exp(-0.5 * sigma^2 * T + sigma .* sqrt(T) .* g) - K, 0);

    optionPrice = B * mean(optionPriceVector);
    optionPrice_COMB = B/2 * mean(optionPriceVector + optionPriceVectorAV);

    optionPricesMC = [optionPricesMC optionPrice];
    optionPricesMC_COMB = [optionPricesMC_COMB optionPrice_COMB];
end 

unbiasedStd = std(optionPricesMC);
unbiasedStdAV = std(optionPricesMC_COMB);

end %function AntitheticVariablesMC