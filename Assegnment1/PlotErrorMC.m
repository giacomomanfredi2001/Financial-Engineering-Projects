function [M,stdEstim]=PlotErrorMC(F0,K,B,T,sigma)
% Plot Error with CRR approach
%
%INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% Given by text, M = 2^m with m = 1:10
%
%OUTPUT
% M:            row-vector of number of simulations
% stdEstim:     row-vector of std deviations of MC simulations
% 
%USES
% function EuropeanOptionPrice(F0, K, B, T, sigma, pricingMode, N, flag)

basisPoint = 1e-4;

Nsim = 100;
M = 2.^(1:20);

stdEstim = [];

for i = 1:length(M)
    optionPricesMC = [];
    
    for j = 1:Nsim
        optionPricesMC = [optionPricesMC EuropeanOptionPrice(F0, K, B, T, sigma, 3, M(i), 1)];
    end 

    UnbiasedSD = std(optionPricesMC);

    stdEstim = [stdEstim UnbiasedSD];
end

figure;
plot(log(M), log(stdEstim), '-*'); hold on;
plot(log(M), log(1./sqrt(M)), '--', 'LineWidth', 2); hold on;
yline(log(basisPoint));

xlabel('Log(M)'); ylabel('Log(stdEstim)');
title('Behaviour of stdEstim'); grid on;
legend('errorCRR', 'log(1/sqrt(M))', 'basisPoint threshold');

p = polyfit(log(M), log(stdEstim), 1);
slope = p(1);
disp(['Estimated Convergence Rate for MC: ', num2str(slope)]);

end %function PlotErrorMC