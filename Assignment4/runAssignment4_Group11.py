# Assignment 4, Group 11
# SW: Marchetto Erica
# Manfredi Giacomo, Maspes Marco, Meschieri Giacomo
import datetime

###################### LIBRARIES #######################

import numpy as np
import pandas as pd
import scipy.io as sio
import FE_Library_A4 as FElib

################### IMPORT OF CSV FILES #########################

dataset = pd.read_csv('EUROSTOXX50_Dataset.csv')

###################  0. Exercise ################################

# Dataset initialization
dataset_0 = dataset.copy()

# Computation of the Log returns
returns = FElib.ReturnComputation(dataset_0, [0, 3, 7, 33, 36], '2015-02-20', '2020-02-20')

# Quantities of interests
alpha = 0.99
notional = 15 * (10**6)
delta = 1

weights = 1/(returns.shape[1]-1) * np.ones(returns.shape[1]-1)

# Computation of the VaR and ES
ES, VaR = FElib.AnalyticalNormalMeasures(alpha, weights.transpose(), notional, delta, returns)

print('E0 - VaR: ', VaR)
print('E0 - ES: ', ES)
print()

###########################################################################

######################### 1. Case Study ###################################

################### POINT A: HS & STATISTICAL BOOTSTRAP ####################

# Historical simulation

# Initialization of the dataset
dataset_1A = dataset.copy()

# Extraction of the stocks and computation of the Log returns
stocks = FElib.StocksExtraction(dataset_1A, [0, 48, 9, 43, 50], '2014-03-20', '2019-03-20')

returns = FElib.ReturnComputation(dataset_1A, [0, 48, 9, 43, 50], '2014-03-20', '2019-03-20')

# Quantities of interest
alpha = 0.95
delta = 1

shares = [25 * 10**3, 20 * 10**3, 20 * 10**3, 10 * 10**3]
values = stocks.iloc[-1, 1:]
portfolioValue = np.dot(shares, values)

weights = (shares*values)/portfolioValue

# Computation of VaR and ES
ES, VaR = FElib.HSMeasurements(returns, alpha, weights, portfolioValue, delta)
print('E1A - VaR Historical Simulation: ', VaR)
print('E1A - ES Historical Simulation: ', ES)

# Bootstrap:

# Quantity of interest
nSim = 200

# Bootstrap to select samples
samples = FElib.bootstrapStatistical(nSim, returns)

# Computation of Var and ES
ES, VaR = FElib.HSMeasurements(samples, alpha, weights, portfolioValue, delta)
print('E1A - VaR Bootstrap: ', VaR)
print('E1A - ES Bootstrap: ', ES)

# Plausibility Check:

VaR = FElib.plausibilityCheck(returns, weights, alpha, portfolioValue, delta)
print('E1A - Plausibility Check: ', VaR)
print()

########################### POINT B: WHS ############################

dataset_1B = dataset.copy()

# Extraction of the stocks and computation of the Log returns
returns = FElib.ReturnComputation(dataset_1B, [0, 3, 5, 12, 13, 19], '2014-03-20', '2019-03-20')

# Quantities of interest
alpha = 0.95
lamda = 0.95
portfolioValue = 1
delta = 1
weights = 1/(returns.shape[1]-1) * np.ones(returns.shape[1]-1)

ES, VaR = FElib.WHSMeasurements(returns, alpha, lamda, weights, portfolioValue, delta)

print('E1B - VaR Weighted Historical Simulation:', VaR)
print('E1B - ES Weighted Historical Simulation: ', ES)

# Plausibility Check:

VaR = FElib.plausibilityCheck(returns, weights, alpha, portfolioValue, delta)
print('E1B - Plausibility Check: ', VaR)
print()

########################### POINT C: PCA ############################

dataset_1C = dataset.copy()
companies = [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

# Extraction of the stocks and computation of the Log returns
stocks = FElib.StocksExtraction(dataset_1C, companies, '2014-03-20', '2019-03-20')
returns = FElib.ReturnComputation(dataset_1C, companies, '2014-03-20', '2019-03-20')

# Quantities of interest
alpha = 0.95
delta = 10
portfolioValue = 1
nMax = 5

weights = 1/(returns.shape[1]-1) * np.ones(returns.shape[1]-1)

businessDays = 256

mean_returns, cov_returns = FElib.MeanCovarianceComputation(returns)

# Initialization of the parameters
VaR = np.zeros(nMax)
ES = np.zeros(nMax)

yearlyMeanReturns = mean_returns * businessDays
yearlyCovariance = cov_returns * businessDays

# Computation of the VaR and ES for PCA
for i in range(1, nMax+1):
    ES[i-1], VaR[i-1] = (
        FElib.PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, delta/businessDays, alpha, i, portfolioValue))

print('E1C - PCA components:', np.arange(1, 6))
print('E1C - VaR PCA:', VaR)
print('E1C - ES PCA: ', ES)

# Plausibility Check:

VaR = FElib.plausibilityCheck(returns, weights, alpha, portfolioValue, delta)
print('E1C - Plausibility Check: ', VaR)
print()

###########################################################################

######################### 2. Exercise ###################################

dataset_2 = dataset.copy()
companies = [0, 13]

# Extraction of the stocks and computation of the Log returns
stocks = FElib.StocksExtraction(dataset_2, companies, '2015-01-16', '2017-01-16')
returns = FElib.ReturnComputation(dataset_2, companies, '2015-01-16', '2017-01-16')

# Extraction of the initial stock
S0 = stocks.iloc[-1, 1]

# Quantities of interest
delta = 10
alpha = 0.95

lamda = 0.95
businessDays = 256

starting_date = datetime.datetime(2017, 1, 16)
finish_date = datetime.datetime(2017, 4, 18)
a = [starting_date, finish_date]

timeToMaturity = FElib.yearfrac(starting_date, finish_date, 3)

# Computation of the number of stocks/calls
notional = 1186680
numberStocks = notional/S0
numberCall = numberStocks

# Parameters call
strike = 25
volatility = 0.154
dividendYield = 0.031
fixedInterestRate = 0.005

# Computation of the Full MonteCarlo
VaR = FElib.FullMonteCarloVaR(returns, numberStocks, numberCall, S0, strike, fixedInterestRate, dividendYield,
                              volatility, lamda, timeToMaturity, delta/businessDays, alpha, businessDays)

print('E2 - VaR Full MonteCarlo evaluation: ', VaR)

# Computation of the Delta Normal Method
VaR = FElib.DeltaApproachVar(returns, numberStocks, numberCall, S0, strike, fixedInterestRate, dividendYield,
                             volatility, lamda, timeToMaturity, delta/businessDays, alpha, businessDays)

print('E2 - VaR Delta Normal evaluation: ', VaR)
print()

#######################################################################

####################### 3. Case Study  ################################

# Uploading of the dates and discounts
# dates = sio.loadmat('dates.mat', squeeze_me=True)
# dates = dates['dates']

discounts = sio.loadmat('discounts.mat', squeeze_me=True)
discounts = discounts['discounts']

# Uploading the dates of the swaps
datesSwaps = sio.loadmat('datesSwaps.mat', squeeze_me=True)
datesSwaps = datesSwaps['datesSwap']

# Quantities of interest
notional = 30 * (10**6)

nYears = 7;
delta = 1

S0 = 1
volatility = 0.2
partCoeff = 0.99
recovery = 0.4

nSimulations = 10**4

# Initialization of the vectors
simulatedStocks = np.zeros(nYears + 1)
simulatedStocks[0] = S0

payoffs = np.zeros(nYears)

# Required values from previous assignments
dates = [datetime.datetime(2008, 2, 19), datetime.datetime(2009, 2, 19),
         datetime.datetime(2010, 2, 19), datetime.datetime(2011, 2, 21),
         datetime.datetime(2012, 2, 20), datetime.datetime(2013, 2, 19),
         datetime.datetime(2014, 2, 19), datetime.datetime(2015, 2, 19)]

# intensities = [0.00482169, 0.00582197, 0.00687360, 0.00867073, 0.00737483, 0.00612470, 0.00854965]
survival_prob = [0.99518991, 0.98939699, 0.98258258, 0.97412284, 0.96696527, 0.96106101, 0.95287929]
discounts = [1, 0.96134578, 0.92689669, 0.89161515, 0.85599116, 0.81988710, 0.78375395, 0.74780268]

# Computation of the NPV Defaultable Approx
exactNPVdef, upperQuantileNPVdef, lowerQuantileNPVdef = (
    FElib.cliquetPrice(np.array(discounts), np.array(survival_prob), dates, nSimulations,
                       nYears, volatility, notional, partCoeff, recovery, 1))

print('E3 - Price Cliquet 7y - Defaultable approx: ', exactNPVdef)
print('E3 - IC Cliquet 7y - Defaultable approx: [', lowerQuantileNPVdef, ',', upperQuantileNPVdef, ']')
print()

# Computation of the NPV Not Defaultable
notDefNPV, upperQuantileNPVNotDef, lowerQuantileNPVNotDef = (
    FElib.cliquetPrice(np.array(discounts), np.array(survival_prob), dates, nSimulations,
                       nYears, volatility, notional, partCoeff, recovery, 3))

print('E3 - Price Cliquet 7y - Not defaultable approx: ', notDefNPV)
print('E3 - IC Cliquet 7y - Not defaultable approx: [', lowerQuantileNPVNotDef, ',', upperQuantileNPVNotDef, ']')
print()