function [datesCDS, survProbs, intensities] = bootstrapCDS_JT(datesDF, datesCDS, spreadsCDS, intensities, survProbs, recovery)
% Bootstrap of the CDS in Jarrow - Turnbull
% 
%INPUT
% datesDF:         dates of the bootstrapped curve
% datesCDS:        dates of the Swap Curve 
% spreadsCDS:      spreads of the Swap Curve
% intensities:     lambda i for each time bucket [ti; ti+1]
% survProbs:       survival probabilities P(t0, ti) for each datesCDS
% recovery:        recovery rate of the obligor
% 
%OUTPUT:
% datesCDS:        dates of the CDS computed
% survProbs:       survival probabilities P(t0, ti) for each datesCDS
% intensities:     lambda i for each time bucket [ti; ti+1]


    % Computation of the interdelta ACT/365 for the survProbs
    delta_ACT_365 = [0; 1; yearfrac(datesDF(1), datesDF(12:12+5), 3)];
    interdelta_365 = delta_ACT_365(2:end) - delta_ACT_365(1:end-1);

    % The cases of the JT approximation imposes lambda = S/(1 - pi)

    intensities_JT = spreadsCDS/(1 - recovery);

    intensities(1) = intensities_JT(1);

    for i = 2 : length(intensities)
        intensities(i) = i*intensities_JT(i) - sum(intensities(1:i-1));
    end

    for i = 2 : length(survProbs)
        survProbs(i) = survProbs(i-1) * exp(-intensities(i-1)*interdelta_365(i-1));
    end

end %function bootstrapCDS_JT