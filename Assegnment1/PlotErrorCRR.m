function [M,errorCRR]=PlotErrorCRR(F0,K,B,T,sigma)
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
% M:            row-vector of number of time steps
% errorCRR:     row-vector of approximation errors
% 
%USES
% function EuropeanOptionPrice(F0, K, B, T, sigma, pricingMode, N, flag)

basisPoint = 1e-4;

M = 2.^(1:10);
errorCRR = [];

for i = 1:length(M)
    errorCRR = [errorCRR abs(EuropeanOptionPrice(F0, K, B, T, sigma, 2, M(i), 1) - ...
        EuropeanOptionPrice(F0, K, B, T, sigma, 1, M(i), 1))];
end

figure;
plot(log(M), log(errorCRR), '-*'); hold on;
plot(log(M), log(1./sqrt(M)), '--', 'LineWidth', 2); hold on;
yline(log(basisPoint));

xlabel('Log(M)'); ylabel('log(ErrorCRR)');
title('Behaviour of errorCRR'); grid on;
legend('errorCRR', 'log(1/M)', 'basisPoint threshold');

p = polyfit(log(M), log(errorCRR), 1);
slope = p(1);
disp(['Estimated Convergence Rate for CRR: ', num2str(slope)]);

end %function PlotErrorCRR