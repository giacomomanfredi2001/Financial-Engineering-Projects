function optionPrice = EuropeanOptionCRR(F0,K,B,T,sigma,N,flag)
%European option price with CRR 
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% N:     number of time steps
% flag:  1 call, -1 put

%% Quantity of interest
dt = T/N;                  % Time step
dx = sigma * sqrt(dt);     % CRR step

u = exp(dx);               % Upward value
q = (1- 1/u)/(u - 1/u);    % Probability upward

r = -log(B)/T;              % Zero rate
B_discounted = exp(-r*dt);  % Discounting factor of each step

%% Computation of the CRR Tree

OptionPriceVector = F0 .* u.^(N:-2:-N);

% Choice of the option
if (flag == 1)
    OptionPriceVector = max(OptionPriceVector - K, 0);
elseif (flag == -1)
    OptionPriceVector = max(K - OptionPriceVector, 0);
end

% Computation of the CRR steps
while (length(OptionPriceVector) > 1)
    OptionPriceVector = B_discounted.*(q*OptionPriceVector(1:end-1) + (1-q)*OptionPriceVector(2:end));
end

% Final value of the option
optionPrice = OptionPriceVector;

end %function EuropeanOptionCRR