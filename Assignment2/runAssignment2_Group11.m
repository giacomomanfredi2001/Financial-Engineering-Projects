% Assignment_2
% Group 11, AY2023-2024
%
%USES
% function readExcelData(filename, formatData)
% function bootstrap(datesSet, ratesSet)
% function zeroRates(dates, discounts)
% function sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts, discounts_DV01)
% function sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)

clear all; close all; clc;
format long;

%% Settings
formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

%% Read market data

% Import data from the Excel dataset

[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

%% Bootstrap

% The function bootstrap compute the Bootstrap for the general dataset
% given from the input Excel file. It gives back the interpolated discount
% factors and related dates. The dates vector includes SettlementDate as first date

[dates, discounts]=bootstrap(datesSet, ratesSet); 

%% Compute Zero Rates

settlementDate = datesSet.settlement;    

% Computation of the zero rates for every discount factor interpolated
zero_rates = zeroRates(dates, discounts)./100;

%% Plot Results

% Conversion in dates to be used in the plot
yf_dates = yearfrac(settlementDate, dates, 3);
settlementDate_conv = datetime(settlementDate, 'ConvertFrom', 'datenum');
correspondingDates = settlementDate_conv + years(yf_dates);

% Deleting the 1st swap for a better visualization
removal_index = length(dates) - length(datesSet.swaps) + 1;
plot_dates = [correspondingDates(1:removal_index-1); correspondingDates(removal_index+1:end)];
plot_discounts = [discounts(1:removal_index-1); discounts(removal_index+1:end)];
plot_zero_rates = [zero_rates(1:removal_index-1); zero_rates(removal_index+1:end)];

% Plot of the Bootstrap results

figure();
yyaxis left
plot(plot_dates, plot_discounts, '-o'); grid on;
xlabel('Dates'); ylabel('B(t0, ti)');

yyaxis right
plot(plot_dates, plot_zero_rates *100, '-*', 'Color', 'red'); grid on;
xlabel('Dates'); ylabel('Zero Rates');
legend('Discount Factors', 'Zero rates');

%%

%% SENSITIVITIES

% Quantities of interests given from the assignment text

fixed_rate = 0.028173;
notional = 1e7;

% Recomputation of the bootstrap for DV01-shifted

% Shifting of the rates of 1bp
ratesSet_modified = ratesSet;
ratesSet_modified.depos = ratesSet_modified.depos + 1e-4;
ratesSet_modified.futures = ratesSet_modified.futures + 1e-4;
ratesSet_modified.swaps = ratesSet_modified.swaps + 1e-4;

% Recomputation of the bootstrap and zero rates after shifting
[dates_modified, discounts_modified]=bootstrap(datesSet, ratesSet_modified);

% Computation of the sensitivities
[DV01, BPV, DV01_z] = sensSwap(settlementDate, datesSet.swaps(1:6), fixed_rate, dates, discounts, discounts_modified)

% Sensitivities for the notional
DV01_total = DV01*notional
BPV_total = BPV*notional
DV01_z_total = DV01_z*notional

%% SENSITIVITIES for IB Fixed Coupon 

% Computation of the Macauly duration
couponPaymentsDates = datesSet.swaps(1:6);
MacD_IB = sensCouponBond(settlementDate, couponPaymentsDates, fixed_rate, dates, discounts)

%%

%% THEORETICAL EXERCISE 1

% Computation of the numerical final formula after the shortcut.

% Quantities of interest
mid_price_S7 = mean(ratesSet.swaps(7, :), 2);

indexB7 = length(discounts) - length(ratesSet.swaps) + 7;
indexB6 = length(discounts) - length(ratesSet.swaps) + 6;
B7 = discounts(indexB7);
B6 = discounts(indexB6);

delta_6_7 = yearfrac(dates(indexB6), dates(indexB7), 6);

% Computation of the price
bond_price = 1 - B7 * (1 + mid_price_S7 * delta_6_7) + B6




