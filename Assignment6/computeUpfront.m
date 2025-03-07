function [prices, spot_vols, upfront] = computeUpfront(flat_vols, strikes, exercised_dates, exercised_discounts, datesSet, yearlyPayments)
% Computation of the upfront of the contract
% 
%INPUT:
% flat_vols:                 flat volatilities from the mkt
% strikes:                   strikes from the market
% exercised_dates:           dates from the bootstrap
% exercised_discounts:       discounts from the bootstrap
% datesSet:                  dates from the market
% yearlyPayments:            payments of the euribor in a year
% 
%OUTPUT:
% prices:                    prices of the caps
% spot_vols:                 euribor spot vols
% upfront:                   upfront of the baseline case
% 
%USES:
% function prices_matrix(flat_volatilities, strikes, exercise_dates, exercise_discounts, settlement_date, yearlyPayments)
% function computation_spot_vols(flat_vols, capsPrices, strikes, exercised_dates, exercised_discounts, settlementDate, number_years, yearlyPayments)
% function couponsComputation (first_coupon, capped_values, base_value, Libor_rates, spot_vols, strikes, yearlyPayments, exercised_dates, exercised_discounts, settlement_date)
% function computationNPV_bank(spol, settlementDate, exercised_dates, exercised_discounts)

    %% Computation of the spot volatilities [0 >> 15y]
    % In this part, we're going to compute the structure of the spot
    % volatilities necessary for our contract, hence we take consideration
    % of each 3 months for the 15 years from the start.

    % Computing caps price from market data
    prices = prices_matrix(flat_vols, strikes, exercised_dates, exercised_discounts, datesSet.settlement, yearlyPayments);
    
    % Computation of the spot volatilities
    
    number_years = [1:10, 12, 15]';
    
    spot_vols = computation_spot_vols(flat_vols,prices,strikes,exercised_dates(1:yearlyPayments * number_years(end)), ...
                    exercised_discounts(1:yearlyPayments * number_years(end)), datesSet.settlement, number_years, yearlyPayments);

    %% Computation of Party B coupons
    % In this section instead we're going to compute the coupons not yet
    % discounted, the coupons are consequences of the caplets/floorlet
    % parity taken from the online papers.

    NPV_invbank = computationNPVB_wo_upfront(datesSet, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments);

    %% Computation of NPV - Upfront
    % Final computation of the upfront comparing the two legs of the
    % contract
    
    % Quantities of interest
    spol = 2 * 1e-2;
    
    % Computation of NPV of Party A (Bank)
    NPV_bank = computationNPV_bank(spol, datesSet.settlement, exercised_dates(1:15*yearlyPayments), exercised_discounts(1:15*yearlyPayments));
    
    % Calibration of the upfront of the Party B (Investment Bank)
    upfront = NPV_bank - NPV_invbank;

end % function computeUpfront