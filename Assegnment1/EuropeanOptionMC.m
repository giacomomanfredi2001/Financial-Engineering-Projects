function optionPrice = EuropeanOptionMC(F0,K,B,T,sigma,N,flag)
%European option price with MC 
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% N:     number of simulations
% flag:  1 call, -1 put

%% Quantity of interests

g = randn(N, 1);      % Normal std distribution

%% Computation of the MC Price

% Choice of the option
if (flag == 1)
    optionPriceVector = max(F0 .* exp(-0.5 * sigma^2 * T - sigma .* sqrt(T) .* g) - K, 0);
elseif (flag == -1)
    optionPriceVector = max(K - F0 .* exp(-0.5 * sigma^2 * T - sigma .* sqrt(T) .* g), 0);
end

% Discounting of the option
optionPrice = B * mean(optionPriceVector);

end %European option price with MC