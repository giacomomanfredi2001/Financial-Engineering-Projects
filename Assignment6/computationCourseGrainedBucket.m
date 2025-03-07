function [sensitivities_certificate, nSwapsHedged] = computationCourseGrainedBucket(datesSet, bucket_sensitivities, matrix_sens_swap, number_contracts)
% Computation of the Course Grained delta sensitivities
% 
%INPUT:
% datesSet:              initial set of dates
% nucket_sensitivities:  sensitivities of the certificate
% matrix_sens_swap:      sensitivities of the swaps
% number_contracts:      # [depos, futures, swaps]
% 
%OUTPUT:
% sensitivities_certificate:          delta_NPV of the certificate
% nSwapsHedged:                       # of swaps to buy/sell for hedging
% 
%USES:
% function computeSingle_DeltaBucket(datesSet_used, ratesSet_used, exercised_dates, flat_vols, strikes, yearlyPayments, upfront)
% function computeNPVswap(settlementDate, datesSet, yf_dates, discounts, fixed_rate)

    %% Quantities of interest

    number_depos = number_contracts(1);
    number_futures = number_contracts(2);

    sensitivities_certificate = zeros(4,1);
    nSwapsHedged = zeros(4,1);

    %% Computation of the 10-15y CG Bucket
    
    % Weights computation
    weights = [zeros(number_depos + number_futures + 9, 1); ...
        interp1([datesSet.swaps(10), datesSet.swaps(15)], [0, 1], datesSet.swaps(10:15));];

    % Computation of Coarse Grained sensitivity certificate
    sensitivities_certificate(4) = weights' * bucket_sensitivities;

    % Computation of the delta NPV swap
    sensitivities_swaps_15y = weights' * matrix_sens_swap(:, 4); 
    
    % Hedging
    nSwapsHedged(4) = -sensitivities_certificate(4)/sensitivities_swaps_15y;

    %% Computation of the 5-10y CG Bucket 
    
    % Weights computation
    weights = [zeros(number_depos + number_futures + 4, 1); ...
        interp1([datesSet.swaps(5), datesSet.swaps(10)], [0, 1], datesSet.swaps(5:10)); ...
        interp1([datesSet.swaps(10), datesSet.swaps(15)], [1, 0], datesSet.swaps(11:15));];

    % Computation of Coarse Grained sensitivity certificate
    sensitivities_certificate(3) = weights' * bucket_sensitivities;

    % Computation of the delta NPV swap
    sensitivities_swaps_15y = weights' * matrix_sens_swap(:, 4);  
    sensitivities_swaps_10y = weights' * matrix_sens_swap(:, 3); 
    
    % Hedging
    nSwapsHedged(3) = -(sensitivities_certificate(3) + nSwapsHedged(4) * sensitivities_swaps_15y) /sensitivities_swaps_10y;
    
    %% Computation of the 2-5y CG Bucket & Hedging

    % Weights computation
    weights = [zeros(number_depos + number_futures + 2, 1); ...
        interp1([datesSet.swaps(2), datesSet.swaps(5)], [0, 1], datesSet.swaps(3:5)); ...
        interp1([datesSet.swaps(5), datesSet.swaps(10)], [1, 0], datesSet.swaps(6:10)); ...
        zeros(5, 1);];

    % Computation of Coarse Grained sensitivity certificate
    sensitivities_certificate(2) = weights' * bucket_sensitivities;

    % Computation of the delta NPV swap
    sensitivities_swaps_15y = weights' * matrix_sens_swap(:, 4);  
    sensitivities_swaps_10y = weights' * matrix_sens_swap(:, 3); 
    sensitivities_swaps_5y = weights' * matrix_sens_swap(:, 2); 
    
    % Hedging
    nSwapsHedged(2) = -(sensitivities_certificate(2) + nSwapsHedged(4) * sensitivities_swaps_15y + ...
        nSwapsHedged(3) * sensitivities_swaps_10y) /sensitivities_swaps_5y;

    %% Computation of the 0-2y CG Bucket & Hedging
    
    weights = [ones(number_depos + number_futures + 2, 1); ...
        interp1([datesSet.swaps(2), datesSet.swaps(5)], [1, 0], datesSet.swaps(3:5)); ...
        zeros(10, 1);];

    % Computation of Coarse Grained sensitivity certificate
    sensitivities_certificate(1) = weights' * bucket_sensitivities;

    % Computation of the delta NPV swap
    sensitivities_swaps_15y = weights' * matrix_sens_swap(:, 4);  
    sensitivities_swaps_10y = weights' * matrix_sens_swap(:, 3); 
    sensitivities_swaps_5y = weights' * matrix_sens_swap(:, 2); 
    sensitivities_swaps_2y = weights' * matrix_sens_swap(:, 1); 

    % Hedging
    nSwapsHedged(1) = -(sensitivities_certificate(1) + nSwapsHedged(4) * sensitivities_swaps_15y + ...
        nSwapsHedged(3) * sensitivities_swaps_10y + nSwapsHedged(2) * sensitivities_swaps_5y) /sensitivities_swaps_2y;
    
end % function computationCourseGrainedBucket