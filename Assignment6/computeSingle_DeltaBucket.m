function [sensitivity, dates, discounts, vec_sens_swap, vec_sens_cap] = computeSingle_DeltaBucket(datesSet_used, ratesSet_used, exercised_dates, ...
    flat_vols, strikes, yearlyPayments, upfront, year_swap, interpolated_spot_vols, strike_CAP)
% Computation of a single sensitivity in order to compute the single
% element of the vector
% 
%INPUT:
% datesSet_used:              set of initial dates
% ratesSet_used.              set of initial discounts
% exercised_dates:            dates of the euribor
% flat_vols:                  matrix of flat volatilities
% strikes:                    vector of strikes
% yearlyPayments:             payments of the euribor each year
% upfront:                    upfront of the certificate
% year_swap:                  year where the bucket is computed onto
% interpolated_spot_vols:       interpolated spot vols for the strike
% strike_CAP:                   interpolated strike for the hedging
% 
%OUTPUT:
% sensitivity:                deltaNPV of the bucket
% dates:                      dates from the bootstrap
% discounts:                  discounts from the bootstrap
% vec_sens_swap:              vector of sensitivities for the swap ith
% 
%USES:
% function bootstrap(datesSet, ratesSet);
% function interpolationDiscounts(dates, discounts, exercised_dates);
% function prices_matrix(flat_volatilities, strikes, exercise_dates, exercise_discounts, settlement_date, yearlyPayments)
% function computation_spot_vols(flat_vols, capsPrices, strikes, exercised_dates, exercised_discounts, settlementDate, number_years, yearlyPayments)
% function computationNPVB_wo_upfront(datesSet, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments)
% function computationNPV_bank(spol, settlementDate, exercised_dates, exercised_discounts)
% function computeNPVswap(settlementDate, dates, discounts, interp_dates, fixed_rate)

    %% Introduction

    vec_sens_swap = zeros(1, 4);

    ratesSet_incomplete = ratesSet_used;

    bps_shift = 1e-4;
    
    %% Conventions

    ACT360 = 2;
    ACT365 = 3;

    %% Bootstrap
    [dates, discounts] = bootstrap(datesSet_used, ratesSet_used);
    
    exercised_discounts = interpolationDiscounts(dates, discounts, exercised_dates);
    
    %% Computing caps price from market data
    prices = prices_matrix(flat_vols, strikes, exercised_dates, exercised_discounts, datesSet_used.settlement, yearlyPayments);
    
    %% Computation of the spot volatilities
    
    number_years = [1:10, 12, 15]';
    
    spot_vols = computation_spot_vols(flat_vols,prices,strikes,exercised_dates(1:yearlyPayments * number_years(end)), ...
                exercised_discounts(1:yearlyPayments * number_years(end)), datesSet_used.settlement, number_years, yearlyPayments);
    
    %% Computation of the NPVs
    NPV_interbank = computationNPVB_wo_upfront(datesSet_used, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments);
    
    spol = 2 * 1e-2;
    
    NPV_bank = computationNPV_bank(spol, datesSet_used.settlement, exercised_dates(1:15*yearlyPayments), exercised_discounts(1:15*yearlyPayments));
    
    %% Computation of the sensitivity
    sensitivity = NPV_bank - (upfront + NPV_interbank);

    %% Computation of the sensitivities of the swaps
    
    vector_years = [2, 5, 10, 15];

    for i = 1:4
        if (vector_years(i) == year_swap)
            vec_sens_swap(i) = computeNPVswap(datesSet_used.settlement, dates, discounts, ....
                    datesSet_used.swaps(1:vector_years(i)), ratesSet_incomplete.swaps(vector_years(i)) - bps_shift);
        else
            vec_sens_swap(i) = computeNPVswap(datesSet_used.settlement, dates, discounts, ....
                    datesSet_used.swaps(1:vector_years(i)), ratesSet_incomplete.swaps(vector_years(i)));
        end
    end

    %% Computation of the caps'prices

    % Useful quantities
    yf_dates = yearfrac(dates(1), exercised_dates, ACT365);
    exercised_discounts = interpolationDiscounts(dates, discounts, exercised_dates);
    
    % Recomputing the Libor rates
    yearfracs = yearfrac(dates(1), exercised_dates, ACT360);
    interdelta = yearfracs(2:end) - yearfracs(1:end-1);
    forward_discounts = exercised_discounts(2:end) ./ exercised_discounts(1:end-1);
    Libor_rates = (1 ./ forward_discounts - 1) ./ interdelta;

    vec_sens_cap = caps_price(exercised_discounts(2 : 5 * yearlyPayments), yf_dates(1 : 5 * yearlyPayments - 1),...
                                   interdelta(1 : 5 * yearlyPayments - 1), Libor_rates(1 : 5 * yearlyPayments - 1),...
                                   strike_CAP, interpolated_spot_vols, 2);

end % function computeSingle_DeltaBucket