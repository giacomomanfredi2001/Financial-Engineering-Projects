function coupons = couponsComputation (first_coupon, capped_values, base_value, Libor_rates, spot_vols, ...
    strikes, yearlyPayments, exercised_dates, exercised_discounts, settlement_date)
% Coupon computation at payment dates
%
%INPUT:
% first_coupon:         initial quarter coupon
% capped_values:        values at which coupons are capped
% base_value:           value added to the Euribor
% Libor_rates:          Libor rates
% spot_vols:            spot volatilities
% strikes:              strikes from the flat volatilities surface
% yearlyPayments:       number of coupon payments per year
% exercised_dates:      dates of coupon exercise
% exercised_discounts:  discounts in exercised dates
% settlement_date:      date of contract settlement
%
%USES:
% function caps_price(discounts_payments, yf_365, interdelta, Libor_rates, strikes, flat_vol, flag)

    %% Time conventions
    
    ACT360 = 2;
    ACT365 = 3;

    %% Sigma spot interpolation
    
    % Computation of 5y strikes
    K = capped_values - ones()*base_value;
    
    % Volatility vector initialization
    sigma_interpolated = zeros(size(spot_vols, 1), 1);
    
    % Interpolation of first 5y
    sigma_interpolated(1:yearlyPayments*5-1) = spline(strikes, spot_vols(1:yearlyPayments*5-1, :), K(1));

    % Interpolation of second 5y
    sigma_interpolated(yearlyPayments*5:yearlyPayments*10-1) = spline(strikes, spot_vols(yearlyPayments*5:yearlyPayments*10-1, :), K(2));
    
    % Interpolation of third 5y
    sigma_interpolated(yearlyPayments*10:yearlyPayments*15-1) = spline(strikes, spot_vols(yearlyPayments*10:yearlyPayments*15-1, :), K(3));
  
   %% Coupons computation
   
   % Quantities of interest
   yearfracs = yearfrac(settlement_date, exercised_dates, ACT360);
   yearfracs365 = yearfrac(settlement_date, exercised_dates, ACT365);
   interdelta = yearfracs(2:end) - yearfracs(1:end-1);
   
   % Coupon vector initialization
   coupons = zeros(size(Libor_rates, 1) + 1, 1);
   
   % First quarter coupon initialization
   coupons(1) = first_coupon * exercised_discounts(1) * yearfracs(1);
   
   % Strikes vectorialization
   strikes_vector = zeros(size(Libor_rates, 1), 1);
   strikes_vector(1:yearlyPayments*5-1) = K(1);
   strikes_vector(yearlyPayments*5:yearlyPayments*10-1) = K(2);
   strikes_vector(yearlyPayments*10:yearlyPayments*15-1) = K(3);
   
   % Caplets computation for the whole vector of volatilities
   caplets = caps_price(exercised_discounts(2:end), yearfracs365(1:end-1), interdelta, Libor_rates, ...
       strikes_vector, sigma_interpolated, 3);
   
   % Computation of floorlet through caplets floorlet parity
   % F = C - [L - K] * delta * B
   floorlets = caplets - (Libor_rates - strikes_vector) .* interdelta .* exercised_discounts(2:end);
   
   % Capped values vectorialization
   capped_vector = zeros(size(Libor_rates, 1), 1);
   capped_vector(1:yearlyPayments*5-1) = capped_values(1);
   capped_vector(yearlyPayments*5:yearlyPayments*10-1) = capped_values(2);
   capped_vector(yearlyPayments*10:yearlyPayments*15-1) = capped_values(3);
   
   % Computation of coupons
   coupons(2:end) = capped_vector .* exercised_discounts(2:end) .* interdelta - floorlets;
   
end % function couponsComputation