import datetime

import numpy as np
import pandas as pd
import math as mt
import scipy.stats as stat


def StocksExtraction(dataset, companies, startingDate, endingDate):
    # Computation of the return matrix
    #
    # INPUT:
    # dataset:         initial dataset with all the companies
    # companies:       list of companies to take in consideration
    # startingDate:    start of period
    # endingDate:      end of period

    # Selection of the companies
    dataset = dataset.iloc[:, companies]

    # Selection of time period
    dataset = dataset[(dataset['Date'] <= endingDate) & (dataset['Date'] >= startingDate)]

    # Fill the NaN values of the dataset
    dataset = dataset.ffill(axis=0)

    return dataset


def ReturnComputation(dataset, companies, startingDate, endingDate):
    # Computation of the return matrix
    #
    # INPUT:
    # dataset:         initial dataset with all the companies
    # companies:       list of companies to take in consideration
    # startingDate:    start of period
    # endingDate:      end of period
    #
    # USES:
    # function StocksExtraction(dataset, companies, startingDate, endingDate)

    # Extraction of the dataset
    dataset = StocksExtraction(dataset, companies, startingDate, endingDate)

    # Computation of Log returns
    returns = dataset.apply(pd.to_numeric, errors='coerce')
    log_returns = np.log(returns.iloc[:, 1:] / returns.iloc[:, 1:].shift(1))
    log_returns = log_returns[1:]
    log_returns.reset_index(drop=True, inplace=True)

    log_returns = pd.concat([dataset.Date[1:].reset_index(drop=True), log_returns], axis=1)

    return log_returns


def MeanCovarianceComputation(returns):
    # Computation of the Mean and returns
    #
    # INPUT:
    # returns:                   daily log returns
    #
    # OUTPUT:
    # mean_returns:              mean of daily log returns
    # cov_returns:               covariance of daily log returns

    mean_returns = returns.iloc[:, 1:].mean()
    cov_returns = pd.DataFrame.cov(returns.iloc[:, 1:])

    return mean_returns, cov_returns


def AnalyticalNormalMeasures(alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay, returns):
    # Computation of ES and VaR simple
    #
    # INPUT:
    # alpha:                               confidence level
    # weights:                             weights of the shares
    # portfolioValue:                      total value of the portfolio
    # riskMeasureTimeIntervalInDay:        lag of the computation
    # returns:                             log returns
    #
    # OUTPUT:
    # VaR:                     value of the VaR
    # ES:                      value of the ES

    # Calculate the mean and variance
    mean_returns = - np.dot(weights.transpose(), returns.iloc[:, 1:].mean())
    cov_returns = pd.DataFrame.cov(returns.iloc[:, 1:])

    std_returns = np.dot(np.dot(weights.transpose(), cov_returns), weights)

    # Calculate the inverse t-student
    nu = 4
    inverted_t_student = stat.t.ppf(alpha, nu)

    # Calculation of the VaR (Value at Risk)
    VaR = (riskMeasureTimeIntervalInDay * mean_returns +
           mt.sqrt(riskMeasureTimeIntervalInDay) * (std_returns**0.5) * inverted_t_student)
    VaR = VaR * portfolioValue

    # Calculation of ES (Expected ShortFall)
    std_ES = (nu + inverted_t_student ** 2) / (nu - 1) * (stat.t.pdf(inverted_t_student, nu)) / (1 - alpha)
    ES = (riskMeasureTimeIntervalInDay * mean_returns +
           mt.sqrt(riskMeasureTimeIntervalInDay) * (std_returns**0.5) * std_ES)
    ES = ES * portfolioValue

    return ES, VaR


def HSMeasurements(returns, alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    # Computation of ES and VaR through the Historical simulation
    #
    # INPUT:
    # returns:                             log returns
    # alpha:                               confidence level
    # weights:                             weights of the shares
    # portfolioValue:                      total value of the portfolio
    # riskMeasureTimeIntervalInDay:        lag of the computation
    #
    # OUTPUT:
    # VaR:                     value of the VaR
    # ES:                      value of the ES

    # Computation of the vector of losses
    Loss = -portfolioValue * np.dot(returns.iloc[:, 1:], weights)

    # Sorting the vector in decreasing order
    orderedLoss = sorted(Loss, reverse=True)

    # VaR computation
    n = len(orderedLoss)
    VaR = np.sqrt(riskMeasureTimeIntervalInDay) * orderedLoss[mt.floor(n*(1-alpha))-1]

    # ES computation
    i = mt.floor(n*(1-alpha))
    selectedLoss = orderedLoss[0: mt.floor(n*(1-alpha))]
    ES = np.sqrt(riskMeasureTimeIntervalInDay) * np.mean(selectedLoss)

    return ES, VaR


def bootstrapStatistical(numberOfSamplesToBootstrap, returns):
    # Computation of the samples for the statistical bootstrap
    #
    # INPUT:
    # returns:                             log returns
    # numberOfSamplesToBootstrap:          required samples

    random_indices = np.random.randint(0, returns.shape[0], size=numberOfSamplesToBootstrap)
    samples = returns.iloc[random_indices, :]

    return samples


def unorderedLossesComputation(returns, weights, lamda, portfolioValue):
    # Computation of the unordered losses
    #
    # INPUT:
    # returns:                 log returns of the companies
    # weights:                 portfolio weights of the companies
    # lamda:                   multiplicative coefficient
    # portfolioValue:          total value of the Portfolio

    # Computation of the vector of losses
    Loss = -portfolioValue * np.dot(returns.iloc[:, 1:], weights)

    # Computation of the time-weights vector
    n = len(Loss)
    c = (1 - lamda) / (1 - lamda ** n)

    exp = np.array(range(n - 1, -1, -1))
    w = c * lamda ** exp

    # Creating the matrix of weights and losses
    unordered_Losses = np.array([[ws, Losses] for ws, Losses in zip(w, Loss)])

    return unordered_Losses


def obtainingLossWHS(unordered_Losses, alpha):
    # Computation of the right index for WHS
    #
    # INPUT:
    # unordered_Losses:            losses not ordered yet
    # alpha:                       confidence level

    # Sorting the matrix in Losses decreasing order
    ordered_indexes = np.transpose(np.argsort(unordered_Losses[:, 1])[::-1])
    orderedLoss = np.take(unordered_Losses, ordered_indexes, axis=0)

    # Searching i* - threshold to choose weights

    i, sum_weights = -1, 0
    while sum_weights <= (1 - alpha):
        i += 1
        sum_weights += orderedLoss[i, 0]

    # Loss computation
    loss = orderedLoss[i, 1]

    return loss, orderedLoss, i


def WHSMeasurements(returns, alpha, lamda, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    # Computation of ES and VaR through the Weighted Historical simulation
    #
    # INPUT:
    # returns:                             log returns
    # alpha:                               confidence level
    # lamda:                               discount parameter
    # weights:                             weights of the shares
    # portfolioValue:                      total value of the portfolio
    # riskMeasureTimeIntervalInDay:        lag of the computation
    #
    # OUTPUT:
    # VaR:                     value of the VaR
    # ES:                      value of the ES

    # Computation of the unordered losses
    unordered_Losses = unorderedLossesComputation(returns, weights, lamda, portfolioValue)

    # VaR computation
    VaR, orderedLoss, index = obtainingLossWHS(unordered_Losses, alpha)
    VaR = np.sqrt(riskMeasureTimeIntervalInDay) * VaR

    # ES computation
    ES = np.dot(orderedLoss[0:index, 0], orderedLoss[0:index, 1])/np.sum(orderedLoss[0:index, 0])
    ES = np.sqrt(riskMeasureTimeIntervalInDay) * ES

    return ES, VaR


def PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, H, alpha,
                      numberOfPrincipalComponents, portfolioValue):
    # Computation of the PCA
    #
    # INPUT:
    # yearlyCovariance:                covariance matrix of the returns
    # yearlyMeanReturns:               mean vector of the returns
    # weights:                         weights of the portfolio
    # H:                               lag of the computation
    # alpha:                           confidence level
    # numberOfPrincipalComponents:     # of components required
    # portfolioValue:                  total value of the portfolio
    #
    # OUTPUT:
    # ES:                      Expected Shortfall of the PCA
    # VaR:                     Value at Risk of the PCA

    # Computation of the eigenvalues (LAMBDA) and eigenvectors(GAMMA)
    eigenvalues, eigenvectors = np.linalg.eig(yearlyCovariance)

    # Order the eigenvalues and eigenvectors
    unorderedMap = np.vstack((eigenvalues.transpose(), eigenvectors))
    indexes = np.argsort(unorderedMap[0, :])[::-1]
    orderedMap = unorderedMap[:, indexes]

    orderedEigenValues = np.diag(orderedMap[0, :])
    orderedEigenVectors = orderedMap[1:, :]

    # Compute the modified vectors
    mu_hat = np.dot(orderedEigenVectors.transpose(), yearlyMeanReturns)
    weights_hat = np.dot(orderedEigenVectors.transpose(), weights)

    # Compute the reduced values
    mu_red = -np.dot(weights_hat[:numberOfPrincipalComponents].transpose(), mu_hat[:numberOfPrincipalComponents])
    sigma2_red = np.dot(np.dot(weights_hat[:numberOfPrincipalComponents].transpose(),
                               orderedEigenValues[:numberOfPrincipalComponents, :numberOfPrincipalComponents]),
                        weights_hat[:numberOfPrincipalComponents])

    # Compute the standard VaR and ES
    std_VaR = stat.norm.ppf(alpha, 0, 1)
    std_ES = stat.norm.pdf(std_VaR)/(1-alpha)

    # Computation of the VaR
    VaR = H * mu_red + mt.sqrt(H*sigma2_red) * std_VaR
    VaR = VaR * portfolioValue

    # Computation of the ES
    ES = H * mu_red + mt.sqrt(H*sigma2_red) * std_ES
    ES = ES * portfolioValue

    return ES, VaR


def plausibilityCheck(returns, portfolioWeights, alpha, portfolioValue, riskMeasureTimeIntervalInDay):
    # Computation of VaR through the Plausibility Check
    #
    # INPUT:
    # returns:                             log returns
    # portfolioWeights:                    weights of the shares
    # alpha:                               confidence level
    # portfolioValue:                      total value of the portfolio
    # riskMeasureTimeIntervalInDay:        lag of the computation
    #
    # OUTPUT:
    # VaR:                     value of the VaR

    # Computation of the correlation matrix
    C = np.corrcoef(returns.iloc[:, 1:], rowvar=False)

    # Computation of the lower and upper bounds
    lower_bound = np.zeros(returns.shape[1] - 1)
    upper_bound = np.zeros(returns.shape[1] - 1)

    for i in range(returns.shape[1] - 1):
        lower_bound[i] = np.quantile(returns.iloc[:, i+1], 1-alpha)
        upper_bound[i] = np.quantile(returns.iloc[:, i+1], alpha)

    # Computation of the signed VaR
    sVaR = 0.5 * portfolioWeights * (np.abs(lower_bound) + np.abs(upper_bound))

    VaR = portfolioValue * np.sqrt(riskMeasureTimeIntervalInDay) * np.sqrt(np.dot(np.dot(sVaR.transpose(), C), sVaR))

    return VaR


def BSPriceCall(F0, B, sigma, K, T):
    # B&S formula for the call
    #
    # INPUT
    # F0:                  initial forward
    # B:                   discount factor
    # sigma:               volatility
    # K:                   strike
    # T:                   time to maturity

    d1 = (np.log(F0 / K) + (0.5 * pow(sigma, 2)) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(F0 / K) - (0.5 * pow(sigma, 2)) * T) / (sigma * np.sqrt(T))
    price = B*(F0 * stat.norm.cdf(d1,0,1) - K * stat.norm.cdf(d2,0,1))

    return price


def FullMonteCarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate,
                      dividend, volatility, underlyingLambda, timeToMaturityInYears, riskMeasureTimeIntervalInYears,
                      alpha, NumberOfDaysPerYears):
    # Computation of the Full Montecarlo Valuation
    #
    # INPUT:
    # logReturns:                            log returns of the asset
    # numberOfShares:                        shares of the asset
    # numberOfCalls:                         number of calls in the portfolio
    # stockPrice:                            initial price of the call
    # strike:                                strike price of the call
    # rate:                                  risk-free rate
    # dividend:                              dividend yield
    # volatility:                            volatility of the model
    # underlyingLambda:                      lambda for the WHS
    # timeToMaturityInYears:                 TTM for the call
    # riskMeasureTimeIntervalInYears:        delta of the VaR
    # alpha:                                 confidence level
    # NumberOfDaysPerYears:                  business days in a year (256)

    # Initialization of the data
    nSimulations = 1000

    total_losses = []
    weights_simulations = []

    # Computation of the unordered losses
    unordered_losses = unorderedLossesComputation(logReturns, [1], underlyingLambda, stockPrice)

    for i in range(nSimulations):
        # Extraction of the losses and log returns related
        indexes = np.random.randint(0, len(logReturns)-1, 10)

        extracted_losses = unordered_losses[indexes]
        extracted_log_returns = logReturns.iloc[indexes, 1]
        extracted_log_returns = np.sum(extracted_log_returns)

        # Weight of the simulation >> correspond to the probability of happening
        weights_simulations.append(np.sum(extracted_losses[:, 0]))

        # Simulation of the stock at (t+delta)
        stockPrices_delta = stockPrice * np.exp(extracted_log_returns)

        # Old Call
        F0 = stockPrice * np.exp((rate - dividend) * timeToMaturityInYears)
        B = np.exp(-rate * timeToMaturityInYears)
        loss_old = BSPriceCall(F0, B, volatility, strike, timeToMaturityInYears)

        # New Call
        F0 = stockPrices_delta * np.exp((rate - dividend) * (timeToMaturityInYears -
                                                             riskMeasureTimeIntervalInYears*NumberOfDaysPerYears/365))
        B = np.exp(-rate * (timeToMaturityInYears - riskMeasureTimeIntervalInYears*NumberOfDaysPerYears/365))

        loss_new = BSPriceCall(F0, B, volatility, strike,
                               timeToMaturityInYears - riskMeasureTimeIntervalInYears*NumberOfDaysPerYears/365)

        loss_call = - numberOfCalls * (loss_new - loss_old)
        loss_stock = - numberOfShares * (stockPrices_delta - stockPrice)

        # Computation of the total loss
        total_losses.append(loss_stock - loss_call)

    # We normalize the weights to 1
    weights_simulations = weights_simulations/np.sum(weights_simulations)

    # Zipping of the computation
    total_losses_weights = np.array([[ws, Losses] for ws, Losses in zip(weights_simulations, total_losses)])

    # Closure of the WHS
    VaR, _, _ = obtainingLossWHS(total_losses_weights, alpha)

    return VaR


def DeltaApproachVar(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate,
                     dividend, volatility, underlyingLambda, timeToMaturityInYears, riskMeasureTimeIntervalInYears,
                     alpha, numberOfDaysPerYear):
    # Computation of the Full Montecarlo Valuation
    #
    # INPUT:
    # logReturns:                            log returns of the asset
    # numberOfShares:                        shares of the asset
    # numberOfCalls:                         number of calls in the portfolio
    # stockPrice:                            initial price of the call
    # strike:                                strike price of the call
    # rate:                                  risk-free rate
    # dividend:                              dividend yield
    # volatility:                            volatility of the model
    # underlyingLambda:                      lambda for the WHS
    # timeToMaturityInYears:                 TTM for the call
    # riskMeasureTimeIntervalInYears:        delta of the VaR
    # alpha:                                 confidence level
    # NumberOfDaysPerYears:                  business days in a year (256)

    # Initialization of the data
    nSimulations = 1000

    total_losses = []
    weights_simulation = []

    # Computation of the unordered losses
    unordered_losses = unorderedLossesComputation(logReturns, [1], underlyingLambda, stockPrice)

    for i in range(nSimulations):
        # Extraction of the losses and log returns related
        indexes = np.random.randint(0, len(logReturns)-1, int(numberOfDaysPerYear*riskMeasureTimeIntervalInYears))

        extracted_losses = unordered_losses[indexes]
        extracted_log_returns = logReturns.iloc[indexes, 1]
        extracted_log_returns = np.sum(extracted_log_returns)

        # Weights of the simulation
        weights_simulation.append(np.sum(extracted_losses[:, 0]))

        # Call
        F0 = stockPrice * np.exp((rate - dividend) * timeToMaturityInYears)
        B = np.exp(-dividend * timeToMaturityInYears)

        d1 = ((np.log(F0 / strike) + (0.5 * volatility * timeToMaturityInYears)) /
              (np.sqrt(volatility * timeToMaturityInYears)))

        sensCall = B * stat.norm.cdf(d1, 0, 1)

        # Computation of the losses
        loss_stock = - stockPrice * extracted_log_returns
        loss_call = - stockPrice * sensCall * extracted_log_returns

        # Computation of the total loss
        total_losses.append(numberOfShares * loss_stock - numberOfCalls * loss_call)

    # Normalization of simulations weights
    weights_simulation = weights_simulation/np.sum(weights_simulation)
    total_losses_weights = np.array([[ws, Losses] for ws, Losses in zip(weights_simulation, total_losses)])

    # Computation of VaR with WHS
    VaR, _, _ = obtainingLossWHS(total_losses_weights, alpha)

    return VaR


def cliquetPrice(discounts, survival_prob, dates, nSimulations, nYears,
                 volatility, notional, partCoeff, recovery, flag):
    # Computation of the cliquet option
    #
    # INPUT:
    # discounts:                     discounts factors B(t0, ti)
    # survival_prob:                 survival probability P(t0, ti)
    # dates:                         reset dates of the option
    # nSimulations:                  total number of simulations required
    # nYears:                        max number of years of the option
    # volatility:                    standard deviation of the option
    # notional:                      total amount of money paid in t0
    # partCoeff:                     parameter of the payoff
    # recovery:                      recovery value in case of default
    # flag:                          1 defaultable, 2 exact defaultable, 3 not defaultable
    #
    # OUTPUT:
    # exactNPVdef:                   price of the option in t0
    # upper_bound:                   IC up of the simulation
    # lower_bound:                   IC down of the simulation

    # Initialization
    NPV = []

    # Computation of the rates [t1 >> t7]
    _, yf365, yf360 = zeroRates(dates, discounts)
    yf365 = np.hstack([0, yf365])
    yf360 = np.hstack([0, yf360])

    interdelta_365 = yf365[1:] - yf365[0:len(yf365) - 1]
    interdelta = yf360[1:] - yf360[0:len(yf360) - 1]

    rates = -np.log(discounts[1:]/discounts[0:len(discounts)-1])/interdelta_365

    # Computation of the year fraction (30/360)
    # yf = yearfrac(pd.to_datetime(dates[0]), pd.to_datetime(dates[1:]), 6)

    match flag:
        case 1:     # Defaultable Bond
            simulated_probs = np.random.uniform(0, 1, nSimulations)

            for i in range(nSimulations):

                # Simulation of the underlying
                underlying = np.ones(nYears+1)
                for j in range(nYears):
                    underlying[j+1] = underlying[j] * np.exp((rates[j] - volatility**2/2)*interdelta[j] +
                                                             volatility * np.sqrt(interdelta[j]) * np.random.normal(0, 1))

                # Computation of the payoff
                payoff = np.maximum(partCoeff*underlying[1:] - underlying[0:len(underlying)-1], 0)

                # Find and extract the indicator of survival
                indicator_obtained = np.where(survival_prob > simulated_probs[i], 1, 0)
                indicator_lost = np.where(survival_prob < simulated_probs[i], 0, 1)

                index = findIndexDefault(indicator_obtained)

                # Compute the NPV
                stock_part = np.dot(payoff, (discounts[1:] * indicator_obtained * survival_prob))
                if index != -1:
                    recovery_part = recovery * np.dot(payoff * indicator_lost,
                                                      (discounts[1:]/discounts[index] * indicator_lost) * discounts[index])
                else:
                    recovery_part = 0

                NPV.append(stock_part + recovery_part)

        case 2: # Exact defaultable
            prova = 1
        case 3: # Not defaultable

            for i in range(nSimulations):

                # Simulation of the underlying
                underlying = np.ones(nYears + 1)

                for j in range(nYears):
                    underlying[j + 1] = underlying[j] * np.exp((rates[j] - volatility ** 2 / 2) * interdelta[j] +
                                                               volatility * np.sqrt(interdelta[j]) * np.random.normal(0, 1))

                # Computation of the payoff
                payoff = np.maximum(partCoeff * underlying[1:] - underlying[0:len(underlying) - 1], 0)

                # Compute the NPV
                stock_part = np.dot(payoff, discounts[1:])
                NPV.append(stock_part)


    # Compute the mean value
    exactNPV = notional * np.mean(NPV)

    # Compute the Confidence Interval
    z_crit = stat.norm.ppf(0.975)
    margin_of_error = notional * z_crit * np.std(NPV) / np.sqrt(len(NPV))
    lower_bound = exactNPV - margin_of_error
    upper_bound = exactNPV + margin_of_error

    return exactNPV, upper_bound, lower_bound


def findIndexDefault(indicator_obtained):
    # Custom function to compute the index of the default time
    #
    # INPUT:
    # indicator_obtained:               indicator function describing all the surv probs greater than default

    # Find the first index that is not 1
    for i in range(len(indicator_obtained)):
        if indicator_obtained[i] != 1:
            return i

    return -1


def yearfrac(date1, date2, basis=0):
    """ Fraction of Year Between Dates.
    This function determines the fraction of a year occurring between two
    dates based on the number days between those dates using a specified
    day count basis. Based on MATLAB's yearfrac function.

    :param date1/date2:     values for dates (in pd.datetime format)
    :param basis:           2 for ACT/360
                            3 for ACT/365
                            6 for 30/360
    :return: year fraction
    """

    if basis == 2:
        return (date2-date1)/360
    elif basis == 3:
        return (date2-date1).days/365
    elif basis == 6:
        d2 = min(date2.day, 30)
        d1 = min(date1.day, 30)
        return (360*(date2.year-date1.year)+30*(date2.month-date1.month)+d2-d1) / 360
    else:
        print("Basis not recognised")
        return None


def zeroRates(dates, discounts):
    # Computation of the Zero rate curve
    #
    # INPUT
    # dates:        vector of all possible dates
    # discounts:    vector of all the discounts related to the date

    # Take vectors as numpy array
    dates = np.array(dates)
    discounts = np.array(discounts)

    # Computation of the year fraction
    yf_dates = yearfrac(pd.to_datetime(dates[0]), pd.to_datetime(dates[1:]), 3)
    yf_dates = np.array(yf_dates)

    # Computation of the year fraction 30/360
    yf_dates_360 = []
    yf_dates_360.append(yearfrac(dates[0], dates[1], 6))
    yf_dates_360.append(yearfrac(dates[0], dates[2], 6))
    yf_dates_360.append(yearfrac(dates[0], dates[3], 6))
    yf_dates_360.append(yearfrac(dates[0], dates[4], 6))
    yf_dates_360.append(yearfrac(dates[0], dates[5], 6))
    yf_dates_360.append(yearfrac(dates[0], dates[6], 6))
    yf_dates_360.append(yearfrac(dates[0], dates[7], 6))

    yf_dates_360 = np.array(yf_dates_360)
    # Computation of the zero rates in percentage units
    zRates = - np.log(discounts[1:])/yf_dates *100

    return zRates, yf_dates, yf_dates_360



