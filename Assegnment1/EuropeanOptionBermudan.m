function optionPrice=EuropeanOptionBermudan(F0,K,B,T,sigma, N)
%European option price with Bermudan Approach
%
%INPUT
% F0:         forward price
% B:          discount factor
% K:          strike
% T:          time-to-maturity
% sigma:      volatility
% N:          number of time steps each month

%% Quantity of interest
dt = T/(N);                 % Time step
dx = sigma * sqrt(dt);      % CRR step

u = exp(dx);                % Upward value
q = (1- 1/u)/(u - 1/u);     % Probability upward

r = -log(B)/T;              % Zero rate
B_discounted = exp(-r*dt);  % Discounting factor of each step

%% Computation of the process

optionPrice = max(F0*u.^(N:-2:-N) - K, 0);

for i=0:N-1 
    optionPrice = B_discounted*(q*optionPrice(1:N-i) + (1-q)*optionPrice(2:N-i+1));

    if (dt*(N-i)==1/12)   
        optionPrice = max(optionPrice, (F0*u.^(N-i-1 : -2 : -N+i+1) - K) );
    elseif ((dt*(N-i)==1/6))
        optionPrice = max(optionPrice, (F0*u.^(N-i-1 : -2 : -N+i+1) - K) ); 
    end
end

end % function EuropeanOptionBermudan