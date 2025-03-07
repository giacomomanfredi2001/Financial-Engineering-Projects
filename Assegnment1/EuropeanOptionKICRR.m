function  optionPrice=EuropeanOptionKICRR(F0,K, KI,B,T,sigma,N) 
%European option price with EU Barrier Up & IN in CRR
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% KI:    barrier
% T:     time-to-maturity
% sigma: volatility
% N:     number of steps

%% Quantity of interest
dt = T/N;                  % Time step
dx = sigma * sqrt(dt);     % CRR step

u = exp(dx);               % Upward value
q = (1- 1/u)/(u - 1/u);    % Probability upward

r = -log(B)/T;              % Zero rate
B_discounted = exp(-r*dt);  % Discounting factor of each step

%% Computation of the CRR Tree

OptionPriceVector = F0 .* u.^(N:-2:-N);
Barrier = OptionPriceVector >= KI;

OptionPriceVector = max(OptionPriceVector.*Barrier - K, 0);

% Computation of the final price
while (length(OptionPriceVector) > 1)
    OptionPriceVector = B_discounted.*(q*OptionPriceVector(1:end-1) + (1-q)*OptionPriceVector(2:end));
end

optionPrice = OptionPriceVector;

end %function EuropeanOptionKICRR