function [upfront_VG, monetary_upfront_VG] = pricing_certificate_VG(datesSet,dates, discounts, ...
    dividend, S0_initial, mktStrikes, mktSurface, M, dz, strike, nSim, spol, notional, yearly_payments)
% Computation of the 2y certificate for the VG
% 
% INPUT:
% datesSet:                       datesSet of the excel file
% dates:                          bootstrapped dates
% discounts:                      bootstrapped discounts
% dividend:                       dividend yield
% S0_initial:                     initial underlying
% mktStrikes:                     vec of strikes from the market
% mktSurface:                     vec of volatilities from the market
% M, dz:                          parameters for the FFT calibration
% strike:                         strike of the triggering option
% nSim:                           number of simulations for the Montecarlo
% spol:                           spread over libor for the swap leg
% notional:                       notional of the certificate
% yearly_payments:                euribor payments in a year
% 
% OUTPUT:
% upfront_VG:              upfront at 2y
% monetary_upfront_VG:     upfront for notional value
% 
% USES:
% function exercised_dates_discounts()
% function interpolationRates()
% function calibration_NIG_VG_parameters()
% function simulation_G()
% function E3_montecarlo()
% function NPV_swap_digital()
% function NPV_structured_coupon()

    alpha_VG = 0; n_years = 2;
    
    % Computation of the dates and discounts
    [payment_dates, reset_dates, payment_discounts, ~] = exercised_dates_discounts(datesSet, dates, discounts, n_years, yearly_payments);
    
    % Interpolation of the required rate at reset date
    r_interpolated = interpolationRates(dates, discounts, reset_dates(1, 1));
    
    % Initial Forward Value for the Quadrature
    B0 = exp(-r_interpolated * reset_dates(1, 2));
    F0 = S0_initial*exp(-dividend * reset_dates(1, 2))/B0;

    % Calibration of the VG parameters
    [sigma_VG, k_VG, eta_VG] = calibration_NIG_VG_parameters(reset_dates(1, 2), B0, F0, mktStrikes, mktSurface, alpha_VG, M, dz, 2);
    
    % Simulation of the MC first year
    G_VG = simulation_G(F0, reset_dates(1, 2), k_VG, nSim, 2);
    [price_call_VG, ~] = E3_montecarlo(strike, sigma_VG, reset_dates(1, 2), k_VG, eta_VG, alpha_VG, B0, F0, nSim, 2, G_VG);
    
    % Check over the probability of trigger
    probability_dig_VG = length(find(price_call_VG == 0))/nSim;
    
    % Creation of a probability vector for better code
    prob_vector_VG = [probability_dig_VG; 1 - probability_dig_VG];
    
    %% Computation of the NPV for both parties
    
    % Computation of Party A: Bank XX
    trigger_years = 1;
    NPV_A_VG = NPV_swap_digital(datesSet.settlement, payment_dates(:, 1), payment_discounts, ...
        prob_vector_VG, spol, yearly_payments, trigger_years);
    
    % Computation of Party B: IB
    
    annual_payment_dates = payment_dates(yearly_payments:yearly_payments:end, 1);
    annual_payment_discounts = payment_discounts(yearly_payments:yearly_payments:end);
    coupon_values = [6*1e-2; 2*1e-2];
    
    NPV_B_VG = NPV_structured_coupon(datesSet.settlement, annual_payment_dates, annual_payment_discounts, ...
        prob_vector_VG, coupon_values);
    
    % Computation of the Upfront
    upfront_VG = NPV_A_VG - NPV_B_VG;
    
    % Computation of the monetary upfront
    monetary_upfront_VG = upfront_VG * notional;

end % function pricing_certificate_VG