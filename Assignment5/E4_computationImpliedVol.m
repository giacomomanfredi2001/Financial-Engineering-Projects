function [model_implied_vol, parameters] = E4_computationImpliedVol(F0, discount, delta, mktStrikes, mktSurface, mktMoneyness, alpha, M, dz, k, r)
% Computation of the implied volatility for the Call surface
% 
% F0:                    initial forward value
% discount:              discount value at maturity
% delta:                 time to maturity
% mktStrikes:            market strikes
% mktSurface:            market volatilities
% mktMoneyness:          market moneyness
% alpha:                 parameter of the model
% M:                     order of FFT steps
% dz:                    time step moneyness FFT
% k:                     volvol of the contract
% r:                     risk free rate
% 
%OUTPUT:
% model_implied_vol:    vector of implied volatilities of the model
% parameters:           calibrated from the market [sigma, k, eta]
% 
%USES:
% function EuropeanOptionClosed(F0,K,B,T,sigma,flag)
% function E3_callPriceLewis(discount, F0, alpha, moneyness, sigma, k, eta, timeToMaturity, flag, M, dz)

    %% Market prices computation
    
    mktPrices = EuropeanOptionClosed(F0, mktStrikes, discount, delta, mktSurface, 1);
    
    %% Model prices
    
    options = optimset('Display', 'Off');
    parameters = lsqnonlin(@(parameters) E3_callPriceLewis(discount, F0, alpha, mktMoneyness, parameters(1), parameters(2), parameters(3), delta, 1, M, dz) - mktPrices, ...
                            [0.001, 0.001, 0.001], [0.001, 0.001, 0.001], [Inf Inf Inf], options);
       
    % Post checking theoretical constraints
    w_signed = (1-alpha)/(k*parameters(1)^2)
    
    if (parameters(3) > - w_signed)
        disp('Theoretical constraints satisfied');
    end
    
    %% Market prices
    
    modelPrices = E3_callPriceLewis(discount, F0, alpha, mktMoneyness, parameters(1), parameters(2), parameters(3), delta, 1, M, dz);
    model_implied_vol = blkimpv(F0, mktStrikes, r, delta, modelPrices);

end % function E4_computationImpliedVol