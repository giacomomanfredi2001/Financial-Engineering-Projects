function [datesCDS, survProbs, intensities] = bootstrapCDS_approx(datesDF, datesCDS, spreadsCDS, discountsCDS, intensities, survProbs, recovery)
% Bootstrap of the CDS in approximated method
% 
%INPUT
% datesDF:         dates of the bootstrapped curve
% datesCDS:        dates of the Swap Curve 
% spreadsCDS:      spreads of the Swap Curve
% discountsCDS:    discounts of the Swap Curve
% intensities:     lambda i for each time bucket [ti; ti+1]
% survProbs:       survival probabilities P(t0, ti) for each datesCDS
% recovery:        recovery rate of the obligor
% 
%OUTPUT:
% datesCDS:        dates of the CDS computed
% survProbs:       survival probabilities P(t0, ti) for each datesCDS
% intensities:     lambda i for each time bucket [ti; ti+1]


    % Computation of the interdelta 30/360
    delta_30_360 = [0; 1; yearfrac(datesDF(1), datesDF(12:12+5), 6)];
    interdelta_30_360 = delta_30_360(2:end) - delta_30_360(1:end-1);

    % Computation of the delta ACT/365 for the intensity step
    delta_ACT_365 = [0; 1; yearfrac(datesDF(1), datesDF(12:12+5), 3)];

    for i = 1: length(datesCDS)
        
        % Computation of default recovery part
        recovery_part = -spreadsCDS(i) * (interdelta_30_360.*discountsCDS)'*survProbs(2:end) + ...       
            (1- recovery)* discountsCDS(1:i-1)'*(survProbs(1:i-1) - survProbs(2:i)) + ...
            (1- recovery)* discountsCDS(i)*survProbs(i);

        % Computation of the multiplicative part
        multiplicative_part = (spreadsCDS(i)*interdelta_30_360(i) + 1 - recovery) *discountsCDS(i);

        % Obtaining of the prob
        survProbs(i+1) = recovery_part/multiplicative_part;

        intensities(i) = - log(survProbs(i+1)/survProbs(i))/(delta_ACT_365(i+1) - delta_ACT_365(i));
    end

end %function bootstrapCDS_approx