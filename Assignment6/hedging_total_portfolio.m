function [nSwapsHedged] = hedging_total_portfolio(datesSet, ratesSet, bps_shift, x_CAP, bucket_sensitivities, matrix_sens_swap, bucket_caps)
% Computation of the Course Grained delta sensitivities after the vega cap
% 
%INPUT:
% datesSet:             initial struct of dates
% ratesSet:             initial struct of rates
% bps_shift:            basis point shifting
% x_CAP:                number of Cap contracts
% bucket_sensitivities: contract sensitivities for each bucket
% matrix_sens_swap:     matrix of swaps sensitivities
% bucket_caps:          sensitivities of Caps
%

    %% Quantities of interest

    % sensitivities_certificate = zeros(4,1);
    nSwapsHedged = zeros(4,1);

    %% Computation of the 10-15y CG Bucket
    
    % Weighs for the coarsed-grained bucket
    ratesSet.swaps(15:end) = ratesSet.swaps(15:end) + bps_shift;
    weights = interp1([datesSet.swaps(10), datesSet.swaps(15)], [0, 1], datesSet.swaps(11:14));
    weights = [0; weights; 1];
   
    % Computation certificate DV01
    DV01_certificate = weights' * bucket_sensitivities(21 : end);

    % Computation swaps 15y DV01
    DV01_swaps15y = weights' * matrix_sens_swap(21 : end,4);
    
    % Hedging
    nSwapsHedged(4) = - DV01_certificate / DV01_swaps15y;

    %% Computation of the 5-10y CG Bucket 
    
    % Weighted shifting
    ratesSet.swaps(10) = ratesSet.swaps(10) + bps_shift;

    weights2 = interp1([datesSet.swaps(10), datesSet.swaps(15)], [1, 0], datesSet.swaps(11:14));
    % ratesSet.swaps(11:14) = ratesSet.swaps(11:14) + weights*bps_shift;

    weights1 = interp1([datesSet.swaps(5), datesSet.swaps(10)], [0, 1], datesSet.swaps(6:9));
    % ratesSet.swaps(6:9) = ratesSet.swaps(6:9) + weights*bps_shift;
    weights = [0;weights1;1;weights2;0];
    
    % Computation DV01 certificate
    DV01_certificate = weights' * bucket_sensitivities(16 : end);

    % Computation DV01 swap 10y
    DV01_swap10y = weights' * matrix_sens_swap(16 : end, 3);

    % Computation DV01 swap 15y
    DV01_swap15y = weights' * matrix_sens_swap(16 : end, 4);

    % Hedging
    nSwapsHedged(3) = -(DV01_certificate + nSwapsHedged(4)*DV01_swap15y)/DV01_swap10y;
    
    %% Computation of the 2-5y CG Bucket & Hedging

    % Weighted shifting
    ratesSet.swaps(5) = ratesSet.swaps(5) + bps_shift;

    weights2 = interp1([datesSet.swaps(5), datesSet.swaps(10)], [1, 0], datesSet.swaps(6:9));

    weights1 = interp1([datesSet.swaps(2), datesSet.swaps(5)], [0, 1], datesSet.swaps(3:4));
    weights = [0; weights1; 1; weights2; 0];
    
    % Computation DV01 certificate
    DV01_certificate = weights' * bucket_sensitivities(13:21);

    % Computation DV01 swap 10y
    DV01_swap10y = weights' * matrix_sens_swap(13:21,3);

    % Computation DV01 swap 15y
    DV01_swap15y = weights' * matrix_sens_swap(13:21,4);
    
    % Computation DV01 swap 5y
    DV01_swap5y = weights' * matrix_sens_swap(13:21,2);

    % Computation DV01 CAP 5Y
    DV01_CAP5y = weights' * bucket_caps(13:21);

    % Hedging
    nSwapsHedged(2) = -(DV01_certificate + nSwapsHedged(4)*DV01_swap15y + ...
        nSwapsHedged(3)*DV01_swap10y + x_CAP*DV01_CAP5y) / DV01_swap5y;

    %% Computation of the 0-2y CG Bucket & Hedging
    
    % Weighted shifting
    
    weights2 = interp1([datesSet.swaps(2), datesSet.swaps(5)], [1, 0], datesSet.swaps(3:4));
    weights1 = ones(13,1);
    weights = [weights1; weights2; 0];

    % Computation DV01 certificate
    DV01_certificate = weights' * bucket_sensitivities(1:16);

    % Computation DV01 swap 10y
    DV01_swap10y = weights' * matrix_sens_swap(1:16,3);

    % Computation DV01 swap 15y
    DV01_swap15y = weights' * matrix_sens_swap(1:16,4);
    
    % Computation DV01 swap 5y
    DV01_swap5y = weights' * matrix_sens_swap(1:16,2);

    % Computation DV01 swap 2y
    DV01_swap2y = weights' * matrix_sens_swap(1:16,1);

    % Computation DV01 CAP 5Y
    DV01_CAP5y = weights' * bucket_caps(1 : 16);

    % Hedging
    nSwapsHedged(1) = -(DV01_certificate + nSwapsHedged(4)*DV01_swap15y + ...
        nSwapsHedged(3)*DV01_swap10y + nSwapsHedged(2)*DV01_swap5y + x_CAP*DV01_CAP5y) / DV01_swap2y;

end % function hedging_total_portfolio