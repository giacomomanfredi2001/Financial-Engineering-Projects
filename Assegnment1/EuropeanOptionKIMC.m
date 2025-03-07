function   optionPrice=EuropeanOptionKIMC(F0,K, KI,B,T,sigma,N)  
%European option price with EU Barrier Up & IN in MC
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% KI:    barrier
% T:     time-to-maturity
% sigma: volatility
% N:     number of simulations

%% Quantity of interest
g = randn(N, 1);      % Normal std distribution

%% Computation of the MC Price

% Computation of the vector of simulations
optionPriceVector = F0 .* exp(-0.5 * sigma^2 * T - sigma .* sqrt(T) .* g);
Barrier = optionPriceVector >= KI;

optionPriceVector = max(optionPriceVector.*Barrier - K, 0);

% Discounting of the final value
optionPrice = B * mean(optionPriceVector);

end %function EuropeanOptionKIMC