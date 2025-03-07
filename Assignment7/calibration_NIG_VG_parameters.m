function [sigma, k, eta] = calibration_NIG_VG_parameters(reset_date, B0, F0, mktStrikes, mktSurface, ...
    alpha, M, dz, flag)
% Calibration of the NIG parameters through the FFT approach
% 
% INPUT:
% dates:           dates bootstrapped
% discounts:       discounts bootstrapped
% reset_date:      date of trigger verification
% dividend:        dividend yield
% mktStrikes:      strikes obtained from the market
% mktSurface:      surface obtained from the market
% alpha:           parameter of the NMVM model    
% M, dz:           parameters of the FFT
% flag:            [1: NIG; 2: VG]
% 
% OUTPUT:
% sigma:           volatility of the process
% k:               volvol of the process
% eta:             skewness of the process
% 
% USES:
% function E3_callPriceLewis(discount, F0, alpha, moneyness, sigma, k, eta, timeToMaturity, flag, M, dz)
    
    %% Choice of the parameters for the NIG
    mktMoneyness = log(F0./mktStrikes);
    mktPrices = EuropeanOptionClosed(F0, mktStrikes, B0, reset_date, mktSurface, 1);

    options = optimset('Display', 'Off');

    parameters = lsqnonlin(@(parameters) E3_callPriceLewis(B0, F0, alpha, mktMoneyness, parameters(1), parameters(2), ...
        parameters(3), reset_date, 1, M, dz, flag) - mktPrices, [0.001, 0.001, 0.001], [0.001, 0.001, 0.001], [Inf Inf Inf], options);
    
    sigma = parameters(1); k = parameters(2); eta = parameters(3);

end % function calibration_NIG_parameters