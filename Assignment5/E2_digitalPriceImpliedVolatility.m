function price = E2_digitalPriceImpliedVolatility(priceBlack, strikesSet, strike, volatilitiesSet, F0, discount, initialDate, finalDate, payoff)
% Computation of digital option price in an implied volatility approach
%
%INPUT
% discount:                    final discount value at maturity
% strike:                      strike price to compute the payoff
% F0:                          initial forward value
% initialDate:                 settlement date 
% finalDate:                   expiry date
% sigma:                       volatility of the contract
% payoff:                      payoff of the digital option

    % Convention
    std_act = 3;
    
    % Interpolation
    % K2 = strikeSet(index), K1 = strike
    
    index = find(strikesSet >= strike, 1);
    sigma = interp1(strikesSet(index-1:index), volatilitiesSet(index-1:index), strike);
    
    % Computation of m
    m = (volatilitiesSet(index) - volatilitiesSet(index-1))/(strikesSet(index) - strikesSet(index-1));

    % Computation of year fraction
    yf = yearfrac(initialDate, finalDate, std_act);

    % Computation of d2
    d1 = log(F0/strike)/(sqrt(yf)*sigma) + 1/2 * sigma * sqrt(yf);
    
    % Computation of Vega with closed formula
    Vega = discount * F0 * sqrt(yf) * normpdf(d1); % per averla in euro
    
    % Computation of price
    price = priceBlack - payoff * m * Vega;

end %function digitalPriceImpliedVolatility