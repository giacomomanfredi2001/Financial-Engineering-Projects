function nCapHedged = computationCourseGrainedVega(datesSet, vega_bucket_certificate, vega_bucket_CAP5y, vega_bucket_CAP15y)
% Computation of the coarse grained vega weighting the flat buckets
% 
%INPUT:
% datesSet:                        initial struct of the dates
% vega_bucket_certificate:         vector containing all the sensitivities vega of the certificate
% vega_bucket_CAP5y:               vector of all cap sensitivities for cap 5y
% vega_bucket_CAP15y:              vector of all cap sensitivities for cap 15y
% 

    %% Quantities of interest
    nCapHedged = zeros(2, 1);

    %% Computing weights
    quoted_years = [1:10 12 15]';
    % Computing weights for 15y
    weights15y = interp1([datesSet.swaps(5), datesSet.swaps(15)], [0, 1], datesSet.swaps(quoted_years(5:end)));
    weights5y = 1 - weights15y;
    weights5y = [1; 1; 1; 1; weights5y];

    %% Hedging 5y - 15y

    % Computation  coarsed-grained DV01 certificate
    DV01_certificate = weights15y' * vega_bucket_certificate(5 : end);

    % Computation coarsed-grained DV01 CAP 15 years
    DV01_CAP_15y = weights15y' * vega_bucket_CAP15y(5 : end);
   
    % Computation position of CAP 15 years 
    nCapHedged(2) = - DV01_certificate / DV01_CAP_15y;     

    %% Hedging 0y - 5y
    
    % Computation  coarsed-grained DV01 certificate
    DV01_certificate = weights5y' * vega_bucket_certificate;

    % Computation coarsed-grained DV01 CAP 5 years
    DV01_CAP_5y = weights5y' * vega_bucket_CAP5y;

    % Computation coarsed-grained DV01 CAP 15 years
    DV01_CAP_15y = weights5y' * vega_bucket_CAP15y;

    % Computation position of CAP 5 years 
    nCapHedged(1) = - (DV01_certificate + nCapHedged(2) * DV01_CAP_15y) / DV01_CAP_5y;

end % function computationCourseGrainedVega