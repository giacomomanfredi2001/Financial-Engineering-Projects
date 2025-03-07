function total_vega = computation_total_vega(flat_vol, strikes, exercised_dates, exercised_discounts, ...
    datesSet, yearlyPayments, upfront)
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
% 
%USES:
% function prices_matrix(flat_volatilities, strikes, exercise_dates, exercise_discounts, settlement_date, yearlyPayments)
% function computation_spot_vols(flat_vols, capsPrices, strikes, exercised_dates, exercised_discounts, settlementDate, number_years, yearlyPayments)
% function computationNPV_bank(spol, settlementDate, exercised_dates, exercised_discounts)
% function computationNPVB_wo_upfront(datesSet, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments)

    %% Quantities of interest

    spol = 2 * 1e-2;
    basis_shift = 1e-4;
    
    %% Bumping of volatilities
    % Modification based on the related case, case total: modification of
    % the whole surface, case bucket: modification only of the required row
    flat_vol_bumped = flat_vol + basis_shift;

    %% Computation of the new caps prices
    caps_price_bumped = prices_matrix(flat_vol_bumped, strikes, exercised_dates, exercised_discounts, datesSet.settlement, yearlyPayments);
    
    %% Computation of the new spot volatilities
    number_years = [1:10 12 15]';
    
    spot_vol_bumped = computation_spot_vols(flat_vol_bumped, caps_price_bumped, strikes, exercised_dates(1 : yearlyPayments * number_years(end)), ...
        exercised_discounts(1 : yearlyPayments * number_years(end)), datesSet.settlement, number_years, yearlyPayments);
    
    %% Computation of the NPVs
    
    % Computation of the NPV_bank
    NPV_bank_bumped = computationNPV_bank(spol, datesSet.settlement, exercised_dates(1 : yearlyPayments * number_years(end)), ...
        exercised_discounts(1 : yearlyPayments * number_years(end)));
    
    % Computation of the NPV_IB
    NPV_IB_bumped = upfront + computationNPVB_wo_upfront(datesSet, exercised_dates, ...
        exercised_discounts, spot_vol_bumped, strikes, yearlyPayments);
    
    % Computation of the total vega as the difference between NPV_ptf_bumped and the normal one
    NPV_ptf_bumped = NPV_bank_bumped - NPV_IB_bumped;
    NPV_ptf = 0;
    
    %% Computation of the vega

    total_vega = NPV_ptf_bumped - NPV_ptf;

end % function computation_vega_sensitivities