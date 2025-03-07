function bermudanPrice=EuropeanOptionBermudanDiv(S0,K,B,T,sigma, N, d)
%European option price with Bermudan Approach with dividends
%
%INPUT
% S0:         initial underlying price
% B:          discount factor
% K:          strike
% T:          time-to-maturity
% sigma:      volatility
% N:          number of time steps each month
% d:          dividends
% 
%USES
% function EuropeanOptionBermudan(F0,K,B,T,sigma, N)
% function EuropeanOptionClosed(F0,K,B,T,sigma,flag)

%% Quantities of interest

bermudanPrice = [];
callPrice = []; 

for i = 1:length(d)           
    F0 = S0*exp(-d(i)*T)/B;           
    bermudanPrice = [bermudanPrice EuropeanOptionBermudan(F0,K,B,T,sigma, N)] ;    
    callPrice = [callPrice EuropeanOptionClosed(F0,K,B,T,sigma, 1)];    
end

% Plotting of the obtained discounted prices
figure;
plot(d, bermudanPrice, '-*'); hold on;
plot(d, callPrice, '-o'); grid on;
title('Dividends Bermudan Option');
xlabel('Dividends'); ylabel('Div values');
legend('Bermudan Prices', 'Call Prices');

end %function EuropeanOptionBermudanDiv