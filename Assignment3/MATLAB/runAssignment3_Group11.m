% Assignment_3
% Group 11, AY2023-2024
% SW: Marchetto
% Manfredi, Maspes, Meschieri
% 
%USES
% function readExcelData( filename, formatData)
% function computationSAW(dates, discounts)
% function bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)
% function pricingFtD(nSim, rho, dates, discounts, datesCDS, spreadsCDS, recovery)

clear all; close all; clc;
format long;

%% Settings

formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

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

%% Convention standards
% Conventions needed to compute the year fractions

std_act_365 = 3;
std_30_360 = 6;

%%

%% 1. Spread Over Asset Swap
% Using the discounting curve of swaps vs euribor 3m, computation of the spread
% ASW of a 3y bond

spreadASW = computationSAW(dates, discounts)

%%

%% 2. CDS Bootstrap: different methods

% Introduction of the data

recovery_ISP = 0.4;                                            % PI value of the ISP
cds_s_ISP = 1e-4*[29; 32; 35; 39; 40; 41];                     % Spread values of the ISP
dates_s_ISP = [datesSet.swaps(1:5); datesSet.swaps(7)];

% Interpolation of the CDS

datesCDS_ISP = datesSet.swaps(1:7);
spreadsCDS_ISP = spline(dates_s_ISP, cds_s_ISP, datesCDS_ISP);

clear dates_s_ISP cds_s_ISP;

%% Bootstrap procedure approximated and exact

[~, survProbs_approx, intensities_approx] = bootstrapCDS(dates, discounts, datesCDS_ISP, spreadsCDS_ISP, 1, recovery_ISP);
[~, survProbs_exact, intensities_exact] = bootstrapCDS(dates, discounts, datesCDS_ISP, spreadsCDS_ISP, 2, recovery_ISP);

%% Bootstrap procedure for JT approximation

[datesCDS, survProbs_JT, intensities_JT] = bootstrapCDS(dates, discounts, datesCDS_ISP, spreadsCDS_ISP, 3, recovery_ISP);

%% Plotting of the values

% Conversion in dates to be used in the plot
yf_dates = yearfrac(datesSet.settlement, datesCDS, std_act_365);
settlementDate_conv = datetime(datesSet.settlement, 'ConvertFrom', 'datenum');
correspondingDates = settlementDate_conv + years(yf_dates);

% Explicit plot of the stairs lambdas
figure();
stairs([settlementDate_conv; correspondingDates], [intensities_approx; intensities_approx(end)], 'LineWidth', 2); hold on;
stairs([settlementDate_conv; correspondingDates], [intensities_exact; intensities_exact(end)], 'LineWidth', 2); hold on;
stairs([settlementDate_conv; correspondingDates], [intensities_JT; intensities_JT(end)], 'LineWidth', 2)
title('Plot of the intensities'); grid on;
xlabel('Dates'); ylabel('Intensity values');
legend('Approximated', 'Exact', 'JT approximation', 'Location', 'southeast');

clear yf_dates settlementDate_conv correspondingDates;

% Difference values between approximation methods

max_difference = max(abs(intensities_exact - intensities_approx))
max_difference_JT = max(abs(intensities_exact - intensities_JT))

%%

%% 3. Price First to Default
% Pricing of a bivariate FtD contract on ISP and UCG through the Li model
% using Gaussian Copulas

% Introduction of the data

recovery_UCG = 0.45;                                        % PI value of the UCG
cds_s_UCG = 1e-4*[34; 39; 45; 46; 47; 47];                  % S values of the UCG
dates_s_UCG = [datesSet.swaps(1:5); datesSet.swaps(7)];  

% Interpolation of the CDS

spreadsCDS_UCG = spline(dates_s_UCG, cds_s_UCG, datesCDS);

clear dates_s_UCG cds_s_UCG;

%% Price of the First to Default

% Quantities of interests                                 
nSim = 40000;                           % Number of simulations
rho = 0.2;                              % Correlation ISP vs UCG

spreadsCDS = [spreadsCDS_ISP spreadsCDS_UCG];
recovery = [recovery_ISP recovery_UCG];

spreadCDS = pricingFtD(nSim, rho, dates, discounts, datesCDS, spreadsCDS, recovery)

%% Plot with multiple rhos
% Plotting of the First to Default spreads with respect to different values
% of the correlation rho

rhos = (0:0.005:0.99);
spreads_rhos = zeros(length(rhos), 1);

for i = 1 : length(rhos)
    spreads_rhos(i) = pricingFtD(nSim, rhos(i), dates, discounts, datesCDS, spreadsCDS, recovery);
end

figure;
plot(rhos, spreads_rhos * 1e4, 'LineWidth', 2); hold on;
plot(rhos(1:20:end), spreads_rhos(1:20:end) * 1e4, 'o', 'Color', 'red', 'LineWidth', 2); grid on;
title('Spreads vs Rhos');
xlabel('Rhos'); ylabel('Spreads FtD');
