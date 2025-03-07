% Assignment 7, Group 10
% 
% SW: Marchetto Erica
% Components: Corti Stefano, Manfredi Giacomo, Maspes Marco
% 
% USES:
% function readExcelData()
% function bootstrap()
% function pricing_certificate_NIG()
% function pricing_certificate_VG()
% function pricing_certificate_Black()
% function exercised_dates_discounts()
% function CouponPut()

clear all; close all; clc;
format long;
warning('off', 'all');

%% Settings & Conventions

ACT360 = 2;
ACT365 = 3;
ACT30360 = 6;

formatData='dd/mm/yyyy'; % Pay attention to your computer settings 

%% Read market data and bootstrapping the IR curve

% Import data from the Excel dataset
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

% The function bootstrap compute the Bootstrap for the general dataset
% given from the input Excel file. It gives back the interpolated discount
% factors and related dates. The dates vector includes SettlementDate as first date
% [1 settlement; 3 depos; 7 futures; 49 swaps]

[dates, discounts] = bootstrap(datesSet, ratesSet); 

% Selection of the €-STOXX-50
struct_STOXX = load('cSelect20230131_B.mat').cSelect;

%%

%% E1: Simulations of the NIG € STOXX 50

%% Quantities of interest

% General 
notional = 100 * 1e6;
yearly_payments = 4; nSim = 1e5;
eps = 10; spol = 1.3 * 1e-2; 

% STOXX 50
S0_initial = struct_STOXX.reference;
dividend = struct_STOXX.dividend;
strike = 3200;

% NIG model
alpha = 1/2;
mktSurface = struct_STOXX.surface;
mktStrikes = struct_STOXX.strikes;

% FFT
M = 19; dz = 0.0049;

%%

%% POINT A: Computation with 1 trigger at 1y
%% POINT D: Computation with 1 ER at 2y, but 2 thresholds

[upfront, monetary_upfront, upfront_3y, monetary_upfront_3y, Ft] = pricing_certificate_NIG(datesSet, dates, discounts, ...
    dividend, S0_initial, mktStrikes, mktSurface, alpha, M, dz, strike, nSim, spol, notional, yearly_payments);

% Printing the values

upfront
monetary_upfront
upfront_3y
monetary_upfront_3y

%% Plot and similar

figure();
histogram(Ft); hold on;
xline(3200, 'Color', 'r', 'LineWidth', 2);
xlabel('Forward value'); ylabel('# simulation');
title('Stoxx50 1y');

%%

%% POINT B: computation of the upfront using the Variance Gamma Model

[upfront_VG, monetary_upfront_VG] = pricing_certificate_VG(datesSet,dates, discounts, ...
    dividend, S0_initial, mktStrikes, mktSurface, M, dz, strike, nSim, spol, notional, yearly_payments)

%% POINT C: 
% It's possible to still compute the price using the variance gamma but
% it's necessary to simulate it again.

%%

%% POINT E: Computation of the upfront using Black using only 2 years

[upfront_black, monetary_upfront_black] = pricing_certificate_Black(datesSet, dates, discounts, ...
    S0_initial, dividend, mktStrikes, mktSurface, strike, eps, spol, yearly_payments, notional)

%% Comparison between models

delta_upfront = upfront - upfront_black
delta_monetary_upfront = monetary_upfront - monetary_upfront_black

%%

%% 2) BERMUDAN SWAPTION PRICING VIA HULL-WHITE

%% Quantities of interest

% Hull-White model parameters
a_HW = 0.11;
sigma_HW = 0.008;

% Strike of the swaption
strike = 5 / 100;

% Bermudan reset dates year fractions
bermudanYearFrac = yearfrac(datesSet.settlement, datesSet.swaps(1:10), ACT365);

% Swap fixed leg year fractions
swapYearFrac = yearfrac(datesSet.settlement, datesSet.swaps(1:10), ACT30360);

% Computing interdelta yearfrac, where the first element is the yearfrac
% from the first year up to the second one
interdeltaSwaps = swapYearFrac(2:end) - swapYearFrac(1:end-1);

% Swaps discounts
yearlyDiscounts = yearlySwapsDiscounts (dates, discounts, datesSet.swaps(1:10));

%% Price computation

% Number of time steps
N = 20000;

% Computing the bermudan swaption price
price = Bermudan_swaption_price (dates, discounts, bermudanYearFrac, interdeltaSwaps, yearlyDiscounts, a_HW, sigma_HW, N, strike)

%% Checking the tree implementation

price_EU = European_swaption_price (dates, discounts, bermudanYearFrac, interdeltaSwaps, yearlyDiscounts, a_HW, sigma_HW, N, strike)

%%

%% JAMSHIDIAN COMPUTATION

% Quantities of interest
strike_CC = 1;
coupon = 5 * 1e-2;
a = 11 * 1e-2; sigma = 0.8 * 1e-2;

yearly_payments = 1; n_years = 10;

% Computation of the Coupon Call

[payment_dates, ~, payment_discounts, ~] = exercised_dates_discounts(datesSet, dates, discounts, n_years, yearly_payments);

% Computation of the forward discounts 
% (cut those up to 2y since we have a non-call 2y)

years = (2:9)';
price_swaptions = zeros(length(years), 1);

for i = 1:length(years)
    price_swaptions(i) = CouponPut(payment_dates, payment_discounts, years(i), yearly_payments, coupon, strike_CC, a, sigma);
end

lower = max(price_swaptions)
price_swaptions
upper = sum(price_swaptions)






















