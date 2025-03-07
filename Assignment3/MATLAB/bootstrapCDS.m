function [datesCDS, survProbs, intensities] = bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)
% Bootstrap of the CDS in different methods
% 
%INPUT
% datesDF:         dates of the bootstrapped curve
% discounts:       discounts of the bootstrapped curve
% datesCDS:        dates of the Swap Curve 
% spreadsCDS:      spreads of the Swap Curve
% flag:            1 (approx), 2 (exact) or 3 (JT)
% recovery:        recovery rate of the obligor
% 
%OUTPUT:
% datesCDS:        dates of the CDS computed
% survProbs:       survival probabilities P(t0, ti) for each datesCDS
% intensities:     lambda i for each time bucket [ti; ti+1]
% 
%USES
% function zeroRates(datesDF, discounts)
% function bootstrapCDS_approx(datesDF, datesCDS, spreadsCDS, discountsCDS, intensities, survProbs, recovery);
% function bootstrapCDS_exact(datesDF, datesCDS, spreadsCDS, discountsCDS, intensities, survProbs, recovery);
% function bootstrapCDS_JT(datesCDS, spreadsCDS, intensities, survProbs, recovery);

    %% Quantities of interest required by all cases

    % Survival probabilities vector:
    % [P(t0, t0); P(t0, t1); .. ; P(t0, t7)] = [1; .... ]
    survProbs = zeros(1 + size(datesCDS, 1), 1);
    survProbs(1) = 1;

    % Intensities lambda (piecewise constat)
    intensities = zeros(size(datesCDS));

    % Year fractions with different conventions
    delta_ACT_365_total = yearfrac(datesDF(1), datesDF(2:end), 3);

    % Discounts for the Swaps' dates 1y >> 7y

    zero_rates = zeroRates(datesDF, discounts)/100;
    zr_1y_swap = interp1(delta_ACT_365_total, zero_rates(2:end), 1, 'linear');
    discount_1y_swap = exp(-zr_1y_swap*1);
    discountsCDS = [discount_1y_swap; discounts(12:12+5)];

    % Clearing of unused variables
    clear zero_rates zr_1y_swap discount_1y_swap;

    %% Different possible cases of simulation

    switch flag
        case 1      % Approximated

            [datesCDS, survProbs, intensities] = bootstrapCDS_approx(datesDF, datesCDS, spreadsCDS, ...
                discountsCDS, intensities, survProbs, recovery);
            
        case 2      % Exact computation

            [datesCDS, survProbs, intensities] = bootstrapCDS_exact(datesDF, datesCDS, spreadsCDS, ...
                discountsCDS, intensities, survProbs, recovery);

        case 3      % JARROW - TURNBULL case

            [datesCDS, survProbs, intensities] = bootstrapCDS_JT(datesDF, datesCDS, spreadsCDS, intensities, survProbs, recovery);

    end

end %function bootstrapCDS
