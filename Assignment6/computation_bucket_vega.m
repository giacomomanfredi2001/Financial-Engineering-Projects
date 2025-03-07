function [vega_bucket_certificate, vega_bucket_CAP5y, vega_bucket_CAP15y] = computation_bucket_vega(flat_vol, strikes, exercised_dates, exercised_discounts, ...
    datesSet, yearlyPayments, upfront, strike_CAP5y, CAP_TTM5y, CAP_initial_price_5y,strike_CAP15y, CAP_TTM15y, CAP_initial_price_15y,index_row)
% Computes the vega sensitivities, both bucket and total
% 
%INPUT:
% flat_vol:                matrix of flat volatilities
% strikes:                 vector of strikes from the mkt
% exercised_dates:         euribor dates
% exercised_discounts:     euribor discounts
% datesSet:                initial structure of the dates
% yearlyPayments:          number of payments per year of euribor
% upfront:                 initial upfront of the floating leg
% index_row:               (ONLY FOR BUCKET) index for the bucket date
% 
%USES:
% function prices_matrix(flat_volatilities, strikes, exercise_dates, exercise_discounts, settlement_date, yearlyPayments)
% function computation_spot_vols(flat_vols, capsPrices, strikes, exercised_dates, exercised_discounts, settlementDate, number_years, yearlyPayments)
% function computationNPV_bank(spol, settlementDate, exercised_dates, exercised_discounts)
% function computationNPVB_wo_upfront(datesSet, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments)

    %% Quantities of interest 
    ACT360 = 2;
    ACT365 = 3;

    yearfracs = yearfrac(datesSet.settlement, exercised_dates, ACT360);
    yearfracs_365 = yearfrac(datesSet.settlement, exercised_dates, ACT365);
    interdelta = yearfracs(2:end) - yearfracs(1:end-1);
    forward_discounts = exercised_discounts(2:end) ./ exercised_discounts(1:end-1);
    Libor_rates = (1 ./ forward_discounts - 1) ./ interdelta;

    spol = 2 * 1e-2;

    %% Bumping of volatilities
    % Modification based on the related case, case total: modification of
    % the whole surface, case bucket: modification only of the required row
    flat_vol_bumped = flat_vol;
    flat_vol_bumped(index_row, :) = flat_vol_bumped(index_row, :) + 1e-4;

    %% Computation of the new caps prices
    caps_price_bumped = prices_matrix(flat_vol_bumped, strikes, exercised_dates, exercised_discounts, datesSet.settlement, yearlyPayments);
    
    %% Computation of the new spot volatilities
    number_years = [1:10 12 15]';
    
    spot_vol_bumped = computation_spot_vols(flat_vol_bumped, caps_price_bumped, strikes, exercised_dates(1 : yearlyPayments * number_years(end)), ...
        exercised_discounts(1 : yearlyPayments * number_years(end)), datesSet.settlement, number_years, yearlyPayments);
    
    %% Computation of the NPVs - certificate
    
    % Computation of the NPV_bank
    NPV_bank_bumped = computationNPV_bank(spol, datesSet.settlement, exercised_dates(1 : yearlyPayments * number_years(end)), ...
        exercised_discounts(1 : yearlyPayments * number_years(end)));
    
    % Computation of the NPV_IB
    NPV_IB_bumped = upfront + computationNPVB_wo_upfront(datesSet, exercised_dates, ...
        exercised_discounts, spot_vol_bumped, strikes, yearlyPayments);
    
    % Computation of the total vega as the difference between NPV_ptf_bumped and the normal one
    NPV_ptf_bumped = NPV_bank_bumped - NPV_IB_bumped;
    NPV_ptf = 0;

    vega_bucket_certificate = NPV_ptf_bumped - NPV_ptf;

    %% Computation of the NPVs - CAP 5Y
    interpolated_spot_vols5y_bumped = spline(strikes, spot_vol_bumped(1 : yearlyPayments * CAP_TTM5y - 1, :),...
                                    strike_CAP5y);

    new_CAP5y = caps_price(exercised_discounts(2 : yearlyPayments * CAP_TTM5y),...
                                   yearfracs_365(1 : yearlyPayments * CAP_TTM5y - 1),...
                                   interdelta(1 : yearlyPayments * CAP_TTM5y - 1),...
                                   Libor_rates(1 : yearlyPayments * CAP_TTM5y - 1),...
                                   strike_CAP5y, interpolated_spot_vols5y_bumped, 2);

    vega_bucket_CAP5y = new_CAP5y - CAP_initial_price_5y;

     %% Computation of the NPVs - CAP 15Y
     interpolated_spot_vols15y_bumped = spline(strikes, spot_vol_bumped(1 : yearlyPayments * CAP_TTM15y - 1, :),...
                                    strike_CAP15y);

     new_CAP15y = caps_price(exercised_discounts(2 : yearlyPayments * CAP_TTM15y),...
                                   yearfracs_365(1 : yearlyPayments * CAP_TTM15y - 1),...
                                   interdelta(1 : yearlyPayments * CAP_TTM15y - 1),...
                                   Libor_rates(1 : yearlyPayments * CAP_TTM15y - 1),...
                                   strike_CAP15y, interpolated_spot_vols15y_bumped, 2);

     vega_bucket_CAP15y = new_CAP15y - CAP_initial_price_15y;
end % function computation_vega_sensitivities