function price = E2_digitalPriceBlack(discount, strike, F0, initialDate, finalDate, sigma, payoff)
% Computation of digital option price according to Black model
%
%INPUT
% discount:                  final discount factor
% strike:                    strike of the model
% F0:                        initial forward value
% initialDate:               settlement date of the contract
% finalDate:                 expiry date of the contract
% sigma:                     volatility of the contract

    % Convention
    std_act = 3;
    
    % Computation of year fraction
    yf = yearfrac(initialDate, finalDate, std_act);

    % Computation of d2
    d2 = log(F0/strike)/(sqrt(yf)*sigma) - 1/2 * sigma * sqrt(yf);
    
    % Computation of digital option price
    price = normcdf(d2) * discount * payoff;

end %function digitalPriceBlack