function [forward_rates, interdelta_fwd] = computationForwardRates(dates, discounts, monitoring_dates)
% Computation of the forward rates
% 
%INPUT:
% dates:                dates from the bootstrap curve
% discounts:            discounts from the bootstrap curve
% monitoring_dates:     dates for the quarterly computation
% 
%OUTPUT:
% forward_rates:        forward rates for the simulation
% interdelta_fwd:       interdelta for the simulation
% 
%USES:
% function zeroRates(dates, discounts)

    %% Conventions
    conv_ACT365 = 3;

    %% Computation of the discounts in MD

    rates = zeroRates(dates, discounts)/100;
    rates_MD = interp1(dates(2:end), rates(2:end), monitoring_dates);

    yf_MD = yearfrac(dates(1), monitoring_dates, conv_ACT365);

    discounts_MD = [exp(-rates_MD .* yf_MD)];

    %% Computation of the forward discounts

    forward_discounts = discounts_MD(2:end)./discounts_MD(1:end-1);

    %% Computation of the forward rates

    interdelta_fwd = yf_MD(2:end) - yf_MD(1:end-1);
    forward_rates = [rates_MD(1); -log(forward_discounts)./interdelta_fwd];

    interdelta_fwd = [yf_MD(1); interdelta_fwd];
    
end % function computationForwardRates

