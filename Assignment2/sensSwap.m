function [DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts, discounts_DV01)
% Computation of the sensitivities DV01, DV01z, BPV
% 
%INPUT
% setDate:                  settlement date
% fixedLegPaymentDates:     payments dates for the fixed leg
% fixedRate:                fixed rate S
% dates:                    dates bootstrapped
% discounts:                discount obtained from the bootstrap
% discounts_DV01:           discount bootstrapped after the shift of 1bp
% 
%OUTPUT:
% DV01:              value of the Dollar Value with shift 1bp
% DV01_z:            value of the Dollar Value with shift of zero rates 1bp
% BPV:               value of the Basis Point Value
% 
%USES:
% function computeNPVswap(settlementDate, datesSet, yf_dates, discounts, fixed_rate)
% function zeroRates(dates, discounts)


    %% DV01

    yf_dates = yearfrac(setDate, dates, 3);

    % Computation of NPV base
    NPV = computeNPVswap(setDate, fixedLegPaymentDates, yf_dates, discounts, fixedRate);
    
    % Computation of NPV Shifted
    
    NPV_shifted = computeNPVswap(setDate, fixedLegPaymentDates, yf_dates, discounts_DV01, fixedRate);
    
    DV01 = abs(NPV_shifted - NPV);
    
    %% DV01 - Z
    
    % Computation of NPV base
    NPV = computeNPVswap(setDate, fixedLegPaymentDates, yf_dates, discounts, fixedRate);
    
    % Computation of NPV shifted
    
    % Shifting of the zero rates
    zero_rates = zeroRates(dates, discounts)./100;
    zero_rates_modified = zero_rates + 1e-4;
    discounts_modified = exp(-yf_dates.*zero_rates_modified);

    NPV_shifted = computeNPVswap(setDate, fixedLegPaymentDates, yf_dates, discounts_modified, fixedRate);

    DV01_z = abs(NPV_shifted - NPV);

    %% BPV

    % Computation of the discount factors
    index = find(dates == fixedLegPaymentDates(end));
    B_swap_vector = discounts(index-length(fixedLegPaymentDates)+1: index);

    % Computation of the inter deltas
    delta_vector = zeros(length(fixedLegPaymentDates)+1, 1);
    delta_vector(1) = setDate;
    delta_vector(2:end) = fixedLegPaymentDates;
    
    delta_yf_swap = yearfrac(delta_vector(1:end-1), delta_vector(2:end), 6);

    % Final computation
    BPV = 1e-4 * delta_yf_swap'*B_swap_vector;

end %function sensSwap