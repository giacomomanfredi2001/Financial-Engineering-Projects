% Assignment_6
% Group 10, AY2023-2024
% SW: Marchetto Erica
% Corti Stefano, Manfredi Giacomo, Maspes Marco
%
%USES
% function readExcelData(filename, formatData)
% function completionSwapsCurve(datesSet, ratesSet, max_year)
% function bootstrap(datesSet, ratesSet)
% function readExcelDataVolatilities(filename)
% funciton interpolationDiscounts(dates, discounts, interp_dates)
% function computeUpfront(flat_vols, strikes, exercised_dates, exercised_discounts, datesSet, yearlyPayments)
% function computeDeltaBuckets(datesSet_incomplete, ratesSet_incomplete, exercised_dates, number_contracts, ...
%     bps_shift, flat_vols, strikes, yearlyPayments, max_year, upfront, interpolated_spot_vols, strike_CAP)
% function computation_total_vega(flat_vol, strikes, exercised_dates, exercised_discounts, ...
%     datesSet, yearlyPayments, upfront)
% function caps_price(discounts_payments, yf_365, interdelta, Libor_rates, strikes, flat_vol, flag)
% function computation_bucket_vega(flat_vol, strikes, exercised_dates, exercised_discounts, ...
%     datesSet, yearlyPayments, upfront, strike_CAP5y, CAP_TTM5y, CAP_initial_price_5y,strike_CAP15y, CAP_TTM15y, CAP_initial_price_15y,index_row)
% function computationCourseGrainedBucket(datesSet, bucket_sensitivities, matrix_sens_swap, number_contracts)
% function vega_hedging(flat_vols,initial_spot_vols,strikes, datesSet, exercised_dates, ...
%     exercised_discounts, yearlyPayments, strike_CAP, CAP_TTM, upfront)
% function hedging_total_portfolio(datesSet, ratesSet,  bps_shift,x_CAP,bucket_sensitivities, matrix_sens_swap, delta_caps)
% function computationCourseGrainedVega(datesSet, vega_bucket_certificate, vega_bucket_CAP5y, vega_bucket_CAP15y)

clear all; close all; clc;
format long;
warning('off', 'all');

%% Settings & Conventions

formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

ACT360 = 2;
ACT365 = 3;

%% Read market data and fixing the tables

% Import data from the Excel dataset
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap_20-2-24', formatData);

% Insert the missing dates for the swap years
initial_dates_swaps = zeros(18, 1);
year_added = [(1:1:12)'; (15:5:30)'; 40; 50];

% Populate the initial swap dates vector

for i = 1:length(initial_dates_swaps)
    initial_dates_swaps(i) = datenum(busdate(datetime(datesSet.settlement, 'ConvertFrom', 'datenum') - ...
        caldays(1) + calyears(year_added(i)), 1, eurCalendar));
end

% Insert the new dates inside the old structure
datesSet.swaps = initial_dates_swaps;

% Compute the mid rates for the available values
ratesSet.swaps = mean(ratesSet.swaps, 2);

datesSet_incomplete = datesSet;
ratesSet_incomplete = ratesSet;

% Completion of the swap curve with all the years in the curve [0 >> 50]
max_year = year_added(end);
[datesSet, ratesSet] = completionSwapsCurve(datesSet, ratesSet, max_year);

%%

%% a) BOOTSTRAP

% The function bootstrap compute the Bootstrap for the general dataset
% given from the input Excel file. It gives back the interpolated discount
% factors and related dates. The dates vector includes SettlementDate as first date
% [1 settlement; 4 depos; 7 futures; 49 swaps]

[dates, discounts] = bootstrap(datesSet, ratesSet); 

%%

%% b) PRICING OF THE UPFRONT

%% Quantities of interest

yearlyPayments = 4;
nYears = 30;

% Computation of the quarterly based dates and discounts in [0 >> 15y]
exercised_dates = zeros(yearlyPayments * nYears , 1);

for i = 1:length(exercised_dates)
    exercised_dates(i) = datenum(busdate(datetime(datesSet.settlement, 'ConvertFrom', 'datenum') ...
        - caldays(1) + calmonths(i*(12/yearlyPayments)), 1, eurCalendar));
end

exercised_discounts = interpolationDiscounts(dates, discounts, exercised_dates);

% Obtainance of the flat volatilities
[flat_vols,strikes] = readExcelDataVolatilities("Caps_vol_20-2-24.xlsx");

% Plotting of the flat volatilities
figure();
[X,Y] = meshgrid([(1:10)'; 12; 15; 20; 25; 30], strikes);
surface_flat_vol = surf(X,Y, flat_vols','FaceAlpha',0.5);
xlabel('Years'); ylabel('Strikes'); zlabel('Volatility surface');
title('Flat volatility surface');

%% Computation of the upfront

[prices, spot_vols, upfront] = computeUpfront(flat_vols, strikes, exercised_dates, exercised_discounts, datesSet, yearlyPayments);

upfront

%% Plots

% Plot of the prices
figure();
[X,Y] = meshgrid(strikes, [(1:10)'; 12; 15; 20; 25; 30]);
surface_prices = surf(X,Y, prices,'FaceAlpha',0.5);
xlabel('Strikes'); ylabel('Years'); zlabel('Cap surface');
title('Cap Prices surface');

% Plot of the spot volatilities
figure();
[X,Y] = meshgrid((0.5:0.25:15)', strikes);
surface_spot_vol = surf(X,Y, spot_vols','FaceAlpha',0.5);
xlabel('Years'); ylabel('Strikes'); zlabel('Cap surface');
title('Spot Volatilities surface');

%%

%% c) DELTA BUCKET DV01
% The delta bucket DV01 modifies the rates in order to get a modification
% for each initial rate (4 depos, 7 futures, 18 swaps) at initial time.

%% Quantities of interest

number_depos = 4;
number_futures = 7;
number_swaps = 15;

bps_shift = 1e-4;

%% Computation of quantities for future hedging

% Computing spot vols for the strike of the CAP
strike_CAP = ratesSet.swaps(5);
interpolated_spot_vols = spline(strikes, spot_vols(1 : yearlyPayments * 5 - 1, :), strike_CAP);

%% Computation of the Delta Bucket Sensitivities
number_contracts = [number_depos, number_futures, number_swaps];

[bucket_sensitivities, matrix_sens_swap, vec_sens_cap] = computeDeltaBuckets(datesSet_incomplete, ratesSet_incomplete, exercised_dates, number_contracts, ...
    bps_shift, flat_vols, strikes, yearlyPayments, max_year, upfront, interpolated_spot_vols, strike_CAP);

bucket_depos = bucket_sensitivities(1:number_depos) *1e4
bucket_futures = bucket_sensitivities(number_depos+1:number_depos+number_futures) *1e4
bucket_swaps = bucket_sensitivities(number_depos+number_futures+1:end) *1e4

figure();
plot([1 7 30 60], bucket_depos); grid on;
xlabel('Days'); ylabel('Delta bucket (bps)');
title('Delta bucket sensitivities depos');

figure();
plot([3 6 9 12 15 18 21], bucket_futures); grid on;
xlabel('Months'); ylabel('Delta bucket (bps)');
title('Delta bucket sensitivities futures');

figure();
plot((1:15), bucket_swaps); grid on;
xlabel('Years'); ylabel('Delta bucket (bps)');
title('Delta bucket sensitivities swaps');

%%

%% d) TOTAL VEGA
% Computation of the Vega modification through the total shifting upward of
% the flat volatilities of 1 bps

% Quantities of interest
notional = 50 * 1e6;

% Computation of the total vega
total_vega = computation_total_vega(flat_vols, strikes, exercised_dates, ...
    exercised_discounts, datesSet, yearlyPayments, upfront)

% Computation of the total vega ptf
vega_ptf = total_vega * notional

%%

%% e) VEGA BUCKET SENSITIVITIES

%% Quantities of interest

strike_CAP5y = ratesSet.swaps(5);
strike_CAP15y = ratesSet.swaps(15);
CAP_TTM5y = 5;
CAP_TTM15y = 15;

%% Computation of initial CAP prices

% Computing Libor rates
yearfracs = yearfrac(datesSet.settlement, exercised_dates, ACT360);
yearfracs_365 = yearfrac(datesSet.settlement, exercised_dates, ACT365);
interdelta = yearfracs(2:end) - yearfracs(1:end-1);
forward_discounts = exercised_discounts(2:end) ./ exercised_discounts(1:end-1);
Libor_rates = (1 ./ forward_discounts - 1) ./ interdelta;

% Computing the spot vols for the strike of the ATM 5Y CAP and ATM 15Y CAP
interpolated_spot_vols5y = spline(strikes, spot_vols(1 : yearlyPayments * CAP_TTM5y - 1, :),...
                                    strike_CAP5y);
interpolated_spot_vols15y = spline(strikes, spot_vols(1 : yearlyPayments * CAP_TTM15y - 1, :),...
                                    strike_CAP15y);
    
% Computing the initial price of the ATM 5Y CAP
CAP_initial_price_5y = caps_price(exercised_discounts(2 : yearlyPayments * CAP_TTM5y),...
                                   yearfracs_365(1 : yearlyPayments * CAP_TTM5y - 1),...
                                   interdelta(1 : yearlyPayments * CAP_TTM5y - 1),...
                                   Libor_rates(1 : yearlyPayments * CAP_TTM5y - 1),...
                                   strike_CAP5y, interpolated_spot_vols5y, 2);
CAP_initial_price_15y = caps_price(exercised_discounts(2 : yearlyPayments * CAP_TTM15y),...
                                   yearfracs_365(1 : yearlyPayments * CAP_TTM15y - 1),...
                                   interdelta(1 : yearlyPayments * CAP_TTM15y - 1),...
                                   Libor_rates(1 : yearlyPayments * CAP_TTM15y - 1),...
                                   strike_CAP15y, interpolated_spot_vols15y, 2);

%% Computation of the vega buckets

% Initialization
vega_bucket_certificate = zeros(size(flat_vols,1) - 3, 1);
vega_bucket_CAP5y = zeros(size(flat_vols,1) - 3, 1);
vega_bucket_CAP15y = zeros(size(flat_vols,1) - 3, 1);

for i = 1 : size(flat_vols,1) - 3
    [vega_bucket_certificate(i),vega_bucket_CAP5y(i),vega_bucket_CAP15y(i)] = computation_bucket_vega(flat_vols, strikes, exercised_dates, exercised_discounts, ...
                                    datesSet, yearlyPayments, upfront, strike_CAP5y, CAP_TTM5y, CAP_initial_price_5y, strike_CAP15y, CAP_TTM15y, CAP_initial_price_15y,i);   
end

comparison_vega_bucket_certificate = sum(vega_bucket_certificate);

% Plot of the vega buckets

figure();
plot([1:10 12 15], vega_bucket_certificate * 1e4); grid on;
xlabel('Years'); ylabel('Vega bucket (bps)');
title('Vega bucket sensitivities');

%%

%% f) COURSE-GRAINED BUCKET (0/2, 2/5, 5/10, 10/15) years
% Computation of the sensitivities in Bps and the number of swaps to hedge
% each course grained bucket

[sensitivities_certificate, nSwapsHedged] = computationCourseGrainedBucket(datesSet, bucket_sensitivities, ...
    matrix_sens_swap, number_contracts);

sensitivities_certificate_notional = sensitivities_certificate * notional

figure();
plot([2 5 10 15], sensitivities_certificate * 1e4); grid on;
xlabel('Years'); ylabel('Delta Coarse bucket (bps)');
title('Delta Coarse bucket sensitivities');

nSwapsHedged_notional = nSwapsHedged * notional

%%

%% g) VEGA AND TOTAL HEDGING

% Quantities of interest
strike_CAP = ratesSet.swaps(5);
CAP_TTM = 5;

number_contracts = [number_depos, number_futures, number_swaps];
bps_shift = 1e-4;

% Hedging the Vega with the ATM 5Y CAP
[CAP_initial_price, x_CAP] = vega_hedging(flat_vols,spot_vols,strikes, datesSet, exercised_dates, exercised_discounts, ...
    yearlyPayments, strike_CAP, CAP_TTM, upfront)

total_cap_position = CAP_initial_price * notional * x_CAP

delta_caps = vec_sens_cap - CAP_initial_price;

% Total hedging with the 2Y, 5Y, 10Y and 15Y swaps
[total_hedging_swaps] = hedging_total_portfolio(datesSet, ratesSet, bps_shift, x_CAP, bucket_sensitivities, matrix_sens_swap, delta_caps);

total_swap_position = notional * total_hedging_swaps

%%

%% h) VEGA BUCKET HEDGE 5y/15y cap

nCapHedged = computationCourseGrainedVega(datesSet, vega_bucket_certificate, vega_bucket_CAP5y, vega_bucket_CAP15y)

total_position_caps_vega_bucket = nCapHedged * notional .* [CAP_initial_price_5y; CAP_initial_price_15y]

