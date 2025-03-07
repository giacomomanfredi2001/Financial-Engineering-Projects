function [bucket_sensitivities, matrix_sens_swaps, vec_sens_caps]  = computeDeltaBuckets(datesSet_incomplete, ratesSet_incomplete, exercised_dates, number_contracts, ...
    bps_shift, flat_vols, strikes, yearlyPayments, max_year, upfront, interpolated_spot_vols, strike_CAP)
% Computation of the DV01 Bucket sensitivities. The bucket are represented
% by each of the original contracts used to generate the bootstrap
% [4 depos, 7 futures, 1y > 12y swap, 15y swap]
% 
%INPUT:
% datesSet_incomplete:          struct of dates to start the bootstrap
% ratesSet_incomplete:          struct of rates to start the bootstrap
% exercised_dates:              dates from 3m to 30 y where the euribor is exercised
% number_contracts:             [#depos, #futures, #swaps]
% bps_shift:                    shift of the ith contract
% flat_vols:                    flat volatilities from the market
% strikes:                      strikes from the market
% yearlyPayments:               payments of the euribor each year
% max_year:                     max year to compute the sensitivities
% upfront:                      upfront of the base case
% interpolated_spot_vols:       interpolated spot vols for the strike
% strike_CAP:                   interpolated strike for the hedging
% 
%OUTPUT:
% bucket_sensitivities:         sensitivities of the certificate
% matrix_sens_swaps:            matrix of the sensitivities of the swaps
% 
%USES:
% function completionSwapsCurve(datesSet, ratesSet, max_year)
% function bootstrap(datesSet, ratesSet)
% function interpolationDiscounts(dates, discounts, interp_dates)
% function prices_matrix(flat_volatilities, strikes, exercise_dates, exercise_discounts, settlement_date, yearlyPayments)
% function computation_spot_vols(flat_vols, capsPrices, strikes, exercised_dates, exercised_discounts, settlementDate, number_years, yearlyPayments)
% function computationNPVB_wo_upfront(datesSet, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments)
% function computationNPV_bank(spol, settlementDate, exercised_dates, exercised_discounts)

    %% Generalization of the incomplete framework
    number_depos = number_contracts(1);
    number_futures = number_contracts(2);
    number_swaps = number_contracts(3);

    bucket_sensitivities = zeros(number_depos + number_futures + number_swaps, 1);
    matrix_sens_swaps = zeros(number_depos+number_futures + number_swaps, 4);
    vec_sens_caps = zeros(number_depos+number_futures + number_swaps, 1);

    % Completion of the curve
    [datesSet_incomplete, ratesSet_incomplete] = completionSwapsCurve(datesSet_incomplete, ratesSet_incomplete, max_year);

    % Initialization
    ratesSet_used = ratesSet_incomplete;
    datesSet_used = datesSet_incomplete;
    
    %% Computation of the modified Depo NPVs
    % We're working with 4 depos so we compute a table for each one and modify
    % its mid rate by adding 1bp.
    
    for i = 1:number_depos
        ratesSet_used.depos(i, :) = ratesSet_incomplete.depos(i, :) + bps_shift;
        
        [bucket_sensitivities(i), ~, ~, matrix_sens_swaps(i, :), ...
            vec_sens_caps(i)] = computeSingle_DeltaBucket(datesSet_used, ratesSet_used, ...
            exercised_dates, flat_vols, strikes, yearlyPayments, upfront, 0, interpolated_spot_vols, strike_CAP);
        
        % Re update of the initial framework
        ratesSet_used = ratesSet_incomplete;
        datesSet_used = datesSet_incomplete;
    end
    
    %% Computation of the modified Futures NPVs
    % We're working with 7 futures so we compute a table for each one and modify
    % its mid rate by adding 1bp.
    
    for i = 1:number_futures
        ratesSet_used.futures(i, :) = ratesSet_incomplete.futures(i, :) + bps_shift;
    
        [bucket_sensitivities(number_depos + i), ~, ~, matrix_sens_swaps(number_depos + i, :), ...
            vec_sens_caps(number_depos + i)] = computeSingle_DeltaBucket(datesSet_used, ratesSet_used, ...
            exercised_dates, flat_vols, strikes, yearlyPayments, upfront, 0, interpolated_spot_vols, strike_CAP);

        % Re update of the initial framework
        ratesSet_used = ratesSet_incomplete;
        datesSet_used = datesSet_incomplete;
    end 
    
    %% Computation of the modified Swaps NPVs
    % We're working with 13 swaps so we compute a table for each one and modify
    % its mid rate by adding 1bp.
    
    for i = 1:number_swaps
        ratesSet_used.swaps(i, :) = ratesSet_incomplete.swaps(i, :) + bps_shift;
    
        [bucket_sensitivities(number_depos + number_futures + i), ~, ~, ...
            matrix_sens_swaps(number_depos + number_futures + i, :), ...
            vec_sens_caps(number_depos + number_futures + i)] = computeSingle_DeltaBucket(datesSet_used, ...
            ratesSet_used, exercised_dates, flat_vols, strikes, yearlyPayments, upfront, i, interpolated_spot_vols, strike_CAP);

        % Re update of the initial framework
        ratesSet_used = ratesSet_incomplete;
        datesSet_used = datesSet_incomplete;
    end 

end % function computeDeltaBuckets