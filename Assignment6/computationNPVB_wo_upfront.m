function NPV = computationNPVB_wo_upfront(datesSet, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments)
% Computation of the NPV of the floating leg
% 
%INPUT:
% datesSet:              initial set of dates
% exercised_dates:       exercised dates of euribor
% exercised_discounts:   exercised discount of euribor
% flat_vols:             matrix of flat volatilities
% strikes:               vector of strikes
% yearlyPayments:        euribor payments in a year
% 
%USES:
% function couponsComputation (first_coupon, capped_values, base_value, Libor_rates, spot_vols, strikes, yearlyPayments, exercised_dates, exercised_discounts, settlement_date)

    %% Convention & Quantities of interest

    ACT360 = 2;

    first_coupon = 3 * 1e-2;
    capped_values = [4.3; 4.6; 5.1] .* 1e-2;
    base_value = 1.1 * 1e-2;
     
    %% Computation of the Libor rates

    yearfracs = yearfrac(datesSet.settlement, exercised_dates(1:15*yearlyPayments), ACT360);
    interdelta = yearfracs(2:end) - yearfracs(1:end-1);
    forward_discounts = exercised_discounts(2:15*yearlyPayments) ./ exercised_discounts(1:15*yearlyPayments-1);
    Libor_rates = (1 ./ forward_discounts - 1) ./ interdelta;
    
    %% Computation of coupons

    coupons = couponsComputation(first_coupon, capped_values, base_value, Libor_rates, spot_vols, ...
        strikes, yearlyPayments, exercised_dates(1:15*yearlyPayments), exercised_discounts(1:15*yearlyPayments), datesSet.settlement);
    
    NPV = sum(coupons);

end % function computationNPVB_wo_upfront