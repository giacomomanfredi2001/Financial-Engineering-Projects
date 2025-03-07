function [upfront_black, monetary_upfront_black] = pricing_certificate_Black(datesSet, dates, discounts, ...
    S0_initial, dividend, mktStrikes, mktSurface, strike, eps, spol, yearly_payments, notional)
% Computation of the 2y certificate for the Black approx
% 
% INPUT:
% datesSet:                       datesSet of the excel file
% dates:                          bootstrapped dates
% discounts:                      bootstrapped discounts
% dividend:                       dividend yield
% S0_initial:                     initial underlying
% mktStrikes:                     vec of strikes from the market
% mktSurface:                     vec of volatilities from the market
% strike:                         strike of the triggering option
% eps:                            noise for the call spread
% spol:                           spread over libor for the swap leg
% notional:                       notional of the certificate
% yearly_payments:                euribor payments in a year
% 
% OUTPUT:
% upfront_black:               upfront of the black in pct
% monetary_upfront_black:      upfront in notional value
% 
% USES:
% function exercised_dates_discounts()
% function interpolationRates()
% function EuropeanOptionClosed()
% function NPV_swap_digital()
% function NPV_structured_coupon()

    % We need to compute the probability of the digital using the black model
    % and consequently applying the same procedure
    
    % Years of the contract
    n_years = 2;
    
    % Computation of the dates and discounts
    [payment_dates, reset_dates, payment_discounts, ~] = exercised_dates_discounts(datesSet, dates, discounts, n_years, yearly_payments);

    % Interpolation of the required rate at reset date
    r_interpolated = interpolationRates(dates, discounts, reset_dates(1, 1));
    
    % Initial Forward Value for the Quadrature
    B0 = exp(-r_interpolated * reset_dates(1, 2));
    F0 = S0_initial*exp(-dividend * reset_dates(1, 2))/B0;

    % Interpolation of the required sigma from the market
    baseline_sigma_blk = interp1(mktStrikes, mktSurface, strike);
    modified_sigma_blk = interp1(mktStrikes, mktSurface, strike + eps);
    
    % Computation of the probability of trigger for the Black model
    baseline_call_blk = EuropeanOptionClosed(F0, strike, B0, reset_dates(1, 2), baseline_sigma_blk, 1);
    modified_call_blk = EuropeanOptionClosed(F0, strike + eps, B0, reset_dates(1, 2), modified_sigma_blk, 1);
    
    probability_dig_blk = 1  - (baseline_call_blk - modified_call_blk)/eps;
    
    % Creation of a vector for better code
    prob_vector_blk = [probability_dig_blk; 1 - probability_dig_blk];
    
    %% Computation of the upfront with Black
    
    % Computation of Party A: Bank XX 
    trigger_years = 1;
    NPV_A_black = NPV_swap_digital(datesSet.settlement, payment_dates(:, 1), payment_discounts, ...
        prob_vector_blk, spol, yearly_payments, trigger_years);
    
    % Computation of Party B: IB
    
    annual_payment_dates = payment_dates(yearly_payments:yearly_payments:end, 1);
    annual_payment_discounts = payment_discounts(yearly_payments:yearly_payments:end);
    coupon_values = [6*1e-2; 2*1e-2];

    NPV_B_black = NPV_structured_coupon(datesSet.settlement, annual_payment_dates, annual_payment_discounts, ...
        prob_vector_blk, coupon_values);
    
    % Computation of the Upfront
    upfront_black = NPV_A_black - NPV_B_black;
    
    % Computation of the monetary upfront
    monetary_upfront_black = upfront_black * notional;

end % function pricing_certificate_Black