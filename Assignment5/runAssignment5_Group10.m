% Assignment 5, Group 10
%
% SW: Marchetto Erica
% Corti Stefano, Manfredi Giacomo, Maspes Marco
% 
%USES:
% function readExcelData(filename, formatData)
% 
% function computationForwardRates(dates, discounts, monitoring_dates)
% function E1_simulationUnderlying(monitoring_rates, dividend, initial_underlying, rho, sigma, interdelta,  nSim);
% function interpolationDiscounts(dates, discounts, interp_disc);
% function E1_computationNPV_bank(spol, protection, settlementDate, exercised_dates, exercised_discounts)
% function E1_computationNPV_investment_bank(upfront, ES, weight, alpha, protection, maturity_discount)
% function IC_floating_leg(ES, weight, protection, maturity_discount)
% 
% function zeroRates(dates, discounts)
% function E2_digitalPriceBlack(discount, strike, F0, initialDate, finalDate, sigma, std_act, payoff)
% function E2_digitalPriceImpliedVolatility(priceBlack, strikesSet, strike, volatilitiesSet, F0, discount, initialDate, finalDate, payoff)
% 
% function E3_callPriceLewis(discount, F0, alpha, moneyness, sigma, k, eta, timeToMaturity, flag, M, dz)
% function E3_montecarlo(x_grid, sigma, TTM, k, eta, alpha, discount, F0, nSim);
% 
% function E4_computationImpliedVol(F0, discount, delta, mktStrikes, mktSurface, mktMoneyness, alpha, M, dz, k, r)

clear all; close all; clc;
format long;

%% Settings

formatData='dd/mm/yyyy'; 
rng(0);

%% Read market data

% Import data from the Excel dataset
[datesSet, ~] = readExcelData('MktData_CurveBootstrap', formatData);

%% Loading of dates and discounts
% The dates are initially a structure, hence we need to load and then
% extract the data from the structure.

% Dates vector is composed by:
% - 1 Settlement
% - 3 IB depos
% - 7 Futures
% - 49 IR Swaps

dates = load("dates.mat");
dates = dates.dates; 

% Discounts are related to the the vector of dates
discounts = load("discounts.mat");
discounts = discounts.discounts;

% Load of the parameters for case Study 2/3/4
dataset = load("cSelect20230131_B.mat");

%% ACT conventions

conv_ACT360 = 2;
conv_ACT365 = 3;
conv_30360 = 6;

%%

%% Case Study 1: Certificate Pricing

%% Quantities of interest

% General parameters
nSim = 100000;

spol = 100 * 1e-4;
upfront = 0.02;
rho = 0.4;

notional = 100 * 1e6;
maturity = 4;                       
protection = 0.95;

% Parameters ENEL
S0_enel = 100;
sigma_enel = 0.161;
div_enel = 0.025;

% Parameters AXA
S0_axa = 200;
sigma_axa = 0.20;
div_axa = 0.027;

%% Monitoring dates and forward rates

monitoring_dates = datesSet.swaps(1:4);

[forward_rates, interdelta_fwd] = computationForwardRates(dates, discounts, monitoring_dates);

%% Simulations of the underlying
% Simulation of the underlying trying a GBM structure

ES = E1_simulationUnderlying(forward_rates, [div_enel; div_axa], [S0_enel; S0_axa], rho, [sigma_enel; sigma_axa], interdelta_fwd,  nSim);

%% Computation of the NPV

% Quantities of interest

% Discount of the 4y Swap
maturity_discount = discounts(14);                 

% Computation of the exercised dates and discounts
exercised_dates = zeros(16, 1);

% Holidays following the Italian calendar
holidays = ['24/03/2008'; '01/05/2008'; '25/12/2008'; '26/12/2008'; '01/01/2009'; ...
            '13/04/2009'; '01/05/2009'; '25/12/2009'; '26/12/2009'; '01/01/2009'; ...
            '05/04/2010'; '01/05/2010'; '25/12/2010'; '26/12/2010'; '01/01/2010'; ...
            '25/04/2011'; '01/05/2011'; '25/12/2011'; '26/12/2011'; '01/01/2011'];

for i = 1:16
    exercised_dates(i) = datenum(busdate(datetime(dates(1), 'ConvertFrom', 'datenum') - caldays(1) + calmonths(i*3), 1, holidays));
end

exercised_discounts = interpolationDiscounts(dates, discounts, exercised_dates);

% Final computation of the partecipation coefficients

value_fixed_leg = E1_computationNPV_bank(spol, protection, dates(1), exercised_dates, exercised_discounts);

options = optimset('Display', 'Off');
partecipation_coefficient = lsqnonlin(@(alpha) nSim*value_fixed_leg - ...
    E1_computationNPV_investment_bank(upfront, ES, 1/2, alpha, protection, maturity_discount), 0, 0, Inf, options)

% Confidence interval
[~, ~, IC] = normfit(IC_floating_leg(ES, 1/2, protection, maturity_discount));
part_coef_IC = (E1_computationNPV_bank(spol, protection, dates(1), exercised_dates, exercised_discounts) - upfront)./flipud(IC)

%%

%% Case Study 2

%% Quantities of interest

notional = 1e7;
alpha = 0.05;
payoff = alpha * notional;

% Extracting dates

initialDate = datesSet.settlement; 
finalDate = datesSet.swaps(1);

% Computation of the year fractions
delta = yearfrac(initialDate, finalDate, conv_ACT365);

% Strike is the reference
strike = dataset.cSelect.reference;

% S0 is equal to the strike since is ATM Spot
S0 = strike;

%% Computation of volatility

% Interpolation of surface values
index = find(dataset.cSelect.strikes >= strike, 1);
volatility = interp1(dataset.cSelect.strikes(index-1:index), dataset.cSelect.surface(index-1:index), strike);

%% Computation of discount

% Computation of zero rates
zRates = zeroRates(dates, discounts);

% Interpolation of zero rates
index = find(dates >= finalDate, 1);
r = interp1(dates(index-1:index), zRates(index-1:index), finalDate)/100;

% Computation of discount
discount = exp(-r*delta);

%% Computation of Forward Price

dividends = dataset.cSelect.dividend;
F0 = S0*exp((r-dividends)*delta);

%% Price with Black Model and with Implied Volatility approach

priceBlack = E2_digitalPriceBlack(discount, strike, F0, initialDate, finalDate, volatility, payoff) 
priceImplVol = E2_digitalPriceImpliedVolatility (priceBlack, dataset.cSelect.strikes, strike, dataset.cSelect.surface, F0, discount, initialDate, finalDate, payoff)

diffPrice = abs(priceBlack-priceImplVol)

%%

%% Case Study 3

%% Quantities of interest

alpha = 1/2;  
nSim = 100000;

sigma = 0.20;
k = 1;
eta = 3;

% F0, discount equal to the Case Study 2 

x_grid = (-0.25: 0.01: 0.25)';

%% Computation of the different approaches

% FFT approach

% Parameters
M = 19;
dz = 0.0048;

price_FFT = E3_callPriceLewis(discount, F0, alpha, x_grid, sigma, k, eta, delta, 1, M, dz);

% Quadrature

price_QUAD = E3_callPriceLewis(discount, F0, alpha, x_grid, sigma, k, eta, delta, 2);

% Montecarlo

simulated_matrix = E3_montecarlo(x_grid, sigma, delta, k, eta, alpha, discount, F0, nSim);
[price_MC, ~, IC] = normfit(simulated_matrix);

%% Plot of the results
% We decided to plot in different figures otherwise it was not possible to
% distinguish the different approaches

figure();
plot(F0.*exp(-x_grid), price_FFT, '--', 'LineWidth', 2); grid on;
xlabel('Strikes'); ylabel('Prices');
title('Call Prices');
legend('FFT approach');

figure();
plot(F0.*exp(-x_grid), price_QUAD, '--', 'LineWidth', 2, 'Color', 'r'); grid on;
xlabel('Strikes'); ylabel('Prices');
title('Call Prices');
legend('Quadrature approach');

figure();
plot(F0.*exp(-x_grid), IC(1, :), '--', 'Color', 'k'); hold on;
plot(F0.*exp(-x_grid), price_MC, '--', 'Color', 'g'); hold on;
plot(F0.*exp(-x_grid), IC(2, :), '--', 'Color', 'k'); grid on;
xlabel('Strikes'); ylabel('Prices');
title('Call Prices');
legend('IC lower', 'MC approach', 'IC upper');

figure();
plot(F0.*exp(-x_grid), price_FFT, '-', 'LineWidth', 3); hold on;
plot(F0.*exp(-x_grid), price_QUAD, '--', 'LineWidth', 2, 'Color', 'r'); hold on;
plot(F0.*exp(-x_grid), price_MC, '-.', 'LineWidth', 2, 'Color', 'g'); grid on;
xlabel('Strikes'); ylabel('Prices');
legend('FFT approach', 'Quadrature approach', 'MC approach');
title('Call Prices - together');



%%  Case with alpha = 2/3

alpha = 2/3;

% FFT approach
price_FFT_mod = E3_callPriceLewis(discount, F0, alpha, x_grid, sigma, k, eta, delta, 1, M, dz);

% Quadrature

price_QUAD_mod = E3_callPriceLewis(discount, F0, alpha, x_grid, sigma, k, eta, delta, 2);

% Plot

figure();
plot(F0.*exp(-x_grid), price_FFT_mod, '--', 'LineWidth', 2); grid on;
xlabel('Strikes'); ylabel('Prices');
title('Call Prices');
legend('FFT approach');

figure();
plot(F0.*exp(-x_grid), price_QUAD_mod, '--', 'LineWidth', 2, 'Color', 'r'); grid on;
xlabel('Strikes'); ylabel('Prices');
title('Call Prices');
legend('Quadrature approach');

%%

%% Case Study 4

%% Quantities of interest

alpha = 1/3;

%% Extraction of data

mktSurface = dataset.cSelect.surface;
mktStrikes = dataset.cSelect.strikes;
mktReference = dataset.cSelect.reference;

mktMoneyness = log(F0./mktStrikes);

%% Computation of the implied volatility

[model_implied_vol, parameters] = E4_computationImpliedVol(F0, discount, delta, mktStrikes, mktSurface, mktMoneyness, alpha, M, dz, k, r);

%% Plot creation

figure();
plot(mktStrikes, model_implied_vol, '--','LineWidth', 2); hold on;
plot(mktStrikes, mktSurface, '--', 'LineWidth', 2); grid on;
title('Implied vs market volatilities'); xlabel('Strikes'); ylabel('Volatilities');
legend('Implied Volatilities', 'Market volatilities');

%% Estimating goodness of calibration

% Computing absolute difference between volatilities
differences = abs(model_implied_vol-mktSurface);

% Computing the maximum relative error (in percentage)
MRE = max(differences./mktSurface)*100






