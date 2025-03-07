function [upfront, monetary_upfront, upfront_3y, monetary_upfront_3y, Ft] = pricing_certificate_NIG(datesSet, dates, discounts, ...
    dividend, S0_initial, mktStrikes, mktSurface, alpha, M, dz, strike, nSim, spol, notional, yearly_payments)
% Computation of the 2y and 3y certificate for the NIG
% 
% INPUT:
% datesSet:                       datesSet of the excel file
% dates:                          bootstrapped dates
% discounts:                      bootstrapped discounts
% dividend:                       dividend yield
% S0_initial:                     initial underlying
% mktStrikes:                     vec of strikes from the market
% mktSurface:                     vec of volatilities from the market
% alpha:                          parameter which describe the NMVM model
% M, dz:                          parameters for the FFT calibration
% strike:                         strike of the triggering option
% nSim:                           number of simulations for the Montecarlo
% spol:                           spread over libor for the swap leg
% notional:                       notional of the certificate
% yearly_payments:                euribor payments in a year
% 
% OUTPUT:
% upfront:              upfront at 2y
% monetary_upfront:     upfront for notional value
% upfront_3y:           upfront at 3y
% monetary_upfront_3y:  upfront for the notional value 3y
% Ft:                   1y simulation of the forward
% 
% USES:
% function exercised_dates_discounts()
% function interpolationRates()
% function calibration_NIG_VG_parameters()
% function simulation_G()
% function E3_montecarlo()
% function NPV_swap_digital()
% function NPV_structured_coupon()
% function NPV_structured_coupon_3y()


    % Years of the contract
    n_years = 2;
    
    % Computation of the dates and discounts
    [payment_dates, reset_dates, payment_discounts, ~] = exercised_dates_discounts(datesSet, dates, discounts, n_years, yearly_payments);
    
    % Interpolation of the required rate at reset date
    r_interpolated = interpolationRates(dates, discounts, reset_dates(1, 1));
    
    % Initial Forward Value for the Quadrature
    B0 = exp(-r_interpolated * reset_dates(1, 2));
    F0 = S0_initial*exp(-dividend * reset_dates(1, 2))/B0;
    
    % Calibration of the NIG parameters
    [sigma, k, eta] = calibration_NIG_VG_parameters(reset_dates(1, 2), B0, F0, mktStrikes, mktSurface, alpha, M, dz, 1);
    
    % Simulation of the MC first year
    G_NIG = simulation_G(F0, reset_dates(1, 2), k, nSim, 1);
    [price_call, Ft] = E3_montecarlo(strike, sigma, reset_dates(1, 2), k, eta, alpha, B0, F0, nSim, 1, G_NIG);
    
    % Check over the probability of trigger
    probability_dig = length(find(price_call == 0))/nSim;
    
    % Creation of a probability vector for better code
    prob_vector = [probability_dig; 1 - probability_dig];
    
    %% Computation of the NPV for both parties
    
    % Computation of Party A: Bank XX
    trigger_years = 1;
    NPV_A = NPV_swap_digital(datesSet.settlement, payment_dates(:, 1), payment_discounts, ...
        prob_vector, spol, yearly_payments, trigger_years);
    
    % Computation of Party B: IB
    
    annual_payment_dates = payment_dates(yearly_payments:yearly_payments:end, 1);
    annual_payment_discounts = payment_discounts(yearly_payments:yearly_payments:end);
    coupon_values = [6*1e-2; 2*1e-2];
    
    NPV_B = NPV_structured_coupon(datesSet.settlement, annual_payment_dates, annual_payment_discounts, ...
        prob_vector, coupon_values);
    
    % Computation of the Upfront
    upfront = NPV_A - NPV_B;
    
    % Computation of the monetary upfront
    monetary_upfront = upfront * notional;

    %%

    %% POINT D: Computation with 2 triggers at 1y and 2y

    % Quantities of interest
    n_years = 3;
    
    % Computation of the dates and discounts
    [payment_dates_3y, reset_dates_3y, payment_discounts_3y, ~] = ...
        exercised_dates_discounts(datesSet, dates, discounts, n_years, yearly_payments);
    
    %% Simulation for the 2nd year
    
    % Interpolation of the required rate at reset date
    r_interpolated_2y = interpolationRates(dates, discounts, reset_dates_3y(2, 1));
    r_interpolated_forward = ((1 + r_interpolated)^reset_dates_3y(1, 2))/((1 + r_interpolated_2y)^reset_dates_3y(2, 2)) - 1;
    B1 = exp(-r_interpolated_forward * (reset_dates_3y(2, 2) - reset_dates_3y(1, 2)));
    
    triggered_indexes = price_call == 0;
    not_triggered_indexes = price_call > 0;

    F1_triggered = Ft(triggered_indexes);
    F1_not_triggered = Ft(not_triggered_indexes);
    
    prob_digital_2y_triggered = zeros(length(F1_triggered), 1);
    prob_digital_2y_not_triggered = zeros(length(F1_not_triggered), 1);
    
    G_NIG = simulation_G(Ft, (reset_dates_3y(2, 2) - reset_dates_3y(1, 2)), k, 1, 1);
    
    for i = 1:length(F1_triggered)
        [price_call, ~] = E3_montecarlo(strike, sigma, (reset_dates_3y(2, 2) - reset_dates_3y(1, 2)), ...
        k, eta, alpha, B1, F1_triggered(i), 1, 1, G_NIG(i));
    
        prob_digital_2y_triggered(i) = (price_call == 0);
    end

    for i = 1:length(F1_not_triggered)
        [price_call, ~] = E3_montecarlo(strike, sigma, (reset_dates_3y(2, 2) - reset_dates_3y(1, 2)), ...
        k, eta, alpha, B1, F1_not_triggered(i), 1, 1, G_NIG(length(F1_triggered)+ i));
    
        prob_digital_2y_not_triggered(i) = (price_call == 0);
    end
    
    prob_digital_2y_triggered = mean(prob_digital_2y_triggered);
    prob_digital_2y_not_triggered = mean(prob_digital_2y_not_triggered);
    
    % Creation of a probability vector for better code
    prob_vector_3y = [probability_dig*prob_digital_2y_triggered; probability_dig*(1 - prob_digital_2y_triggered); ...
        (1 - probability_dig)*prob_digital_2y_not_triggered; (1 - probability_dig)*(1 - prob_digital_2y_not_triggered) ];

    %% Computation of the NPV for both parties
    
    trigger_years_3y = [2; 2; 2];
    NPV_A_3y = NPV_swap_digital(datesSet.settlement, payment_dates_3y(:, 1), payment_discounts_3y, ...
        prob_vector_3y, spol, yearly_payments, trigger_years_3y);
    
    % Computation of Party B: IB
    
    annual_payment_dates_3y = payment_dates_3y(yearly_payments:yearly_payments:end, 1);
    annual_payment_discounts_3y = payment_discounts_3y(yearly_payments:yearly_payments:end);
    coupon_values_3y = [6*1e-2; 2*1e-2];
    
    NPV_B_3y = NPV_structured_coupon_3y(datesSet.settlement, annual_payment_dates_3y, annual_payment_discounts_3y, ...
        prob_vector_3y, coupon_values_3y);
    
    % Computation of the Upfront
    upfront_3y = NPV_A_3y - NPV_B_3y;
    
    % Computation of the monetary upfront
    monetary_upfront_3y = upfront_3y * notional;

end % function pricing_certificate_NIG