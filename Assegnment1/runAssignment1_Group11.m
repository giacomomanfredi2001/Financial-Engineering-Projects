% Assignment_1
% Group 11, AA2023-2024
% 
%USES
% function EuropeanOptionPrice(F0,K,B,T,sigma,pricingMode,N,flag)
% function PlotErrorCRR(F0,K,B,T,sigma)
% function PlotErrorMC(F0,K,B,T,sigma)
% function findBest(M, error, threshold)
% function EuropeanOptionKICRR(F0,K, KI,B,T,sigma,N) 
% function EuropeanOptionKIMC(F0,K, KI,B,T,sigma,N)  
% function EuropeanOptionClosed(F0,K,B,T,sigma,flag)
% function VegaKI(F0,K,KI,B,T,sigma,N,flagNum)
% function AntitheticVariablesMC(F0, K, B, T, sigma, N)
% function EuropeanOptionBermudan(F0,K,B,T,sigma, N)
% function EuropeanOptionBermudanDiv(S0,K,B,T,sigma, N, d)

%% Pricing parameters
S0=1;
K=1;
r=0.03;
TTM=1/4; 
sigma=0.22;
flag=1; % flag:  1 call, -1 put
d=0.06;

N = 1e6; % Number of contracts

%% Quantity of interest
B=exp(-r*TTM); % Discount

%% Pricing with Garman & Kohlhagen
F0=S0*exp(-d*TTM)/B;     % Forward in G&C Model F0 = S0*exp((r-d)*TTM)
M=100;                   % M = simulations for MC, steps for CRR;

OptionPrices = [];       % Vector of prices [BLACK CRR MC]

for pricingMode = 1:3
    OptionPrices = [OptionPrices ...
        EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag)];
end

TotalOptionPrices = 1e6*OptionPrices;  % Prices of the total amount of contracts

%% Errors Rescaling 

basisPoint = 1e-4;

% plot Errors for CRR varing number of steps
% Note: both functions plot also the Errors of interest as side-effect 
[nCRR,errCRR]=PlotErrorCRR(F0,K,B,TTM,sigma);
best_M_CRR = findBest(nCRR, errCRR, basisPoint)
best_M_CRR_error = errCRR(find(nCRR == best_M_CRR))

% plot Errors for MC varing number of simulations N 
[nMC,stdEstim]=PlotErrorMC(F0,K,B,TTM,sigma); 
best_M_MC = findBest(nMC, stdEstim, basisPoint)
best_M_MC_error = stdEstim(find(nMC == best_M_MC))

%% KI Option

barrier = 1.3;    % Value of UP&IN barrier threshold

% Computation of the numerical simulations
priceKICRR_numerical = EuropeanOptionKICRR(F0,K, barrier,B,TTM,sigma,best_M_CRR) 
priceKIMC_numerical = EuropeanOptionKIMC(F0,K, barrier,B,TTM,sigma,best_M_MC) 

% Closed formula
closedCall = EuropeanOptionClosed(F0, barrier, B, TTM, sigma, 1);
d2Barrier = (log(F0/barrier) - TTM*sigma^2/2)/(sigma*sqrt(TTM));
closedDigital = B*(barrier-K)*normcdf(d2Barrier);

priceClosed = closedCall + closedDigital

%% KI Option Vega

% Vector of the S0 changing
S0_vector = (0.7:0.01:1.5);
F0_vector = F0.*S0_vector;

% Computation of the vega
M = [3000, best_M_MC, 0];
vecString = ["CRR", "MC", "Closed Formula"];

for i = 1:3
    vega = [];
    
    for j = 1:length(F0_vector)
        vega = [vega VegaKI(F0_vector(j),K,barrier,B,TTM,sigma,M(i),i)];
    end
    
    figure;
    plot(S0_vector, vega); hold on; grid on;
    title(['Vega ', vecString(i)]);
    xlabel('S0'); ylabel('Vega values');
end

%% Antithetic Variables

[unbiasedStd, unbiasedStdAV] = AntitheticVariablesMC(F0, K, B, TTM, sigma, best_M_MC)

%% Bermudan Options

months = 3;         % Months for the pricing
interstep = 6;      % Interstep for each months
N = months * interstep;

% Computation of the Bermudan Price
bermudanPrice = EuropeanOptionBermudan(F0,K,B,TTM,sigma, N)

%% Bermudan Options with dividends

% Quantities of interest
N= 51;  
d = 0:0.005:0.06;      

bermudan_dividends = EuropeanOptionBermudanDiv(S0,K,B,TTM,sigma, N, d)
