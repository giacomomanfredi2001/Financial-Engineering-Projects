import pandas as pd
import numpy as np
import scipy as sc
import keras as k
import tensorflow as tf
import matplotlib.pyplot as plt
import statsmodels.api as sm
from FE_Library import yearfrac
from FE_Library import minimizationFunction

# ------------------------------ Load Time Series ------------------------------------
quotes = pd.read_csv("DatasetPythonAss3.csv")
quotes["Date"] = pd.to_datetime(quotes["Date"], format="%d/%m/%Y")
quotes = quotes.set_index("Date")
AAPL = quotes["AAPL"]
SPX = quotes["SPX"]

# -------------------------------  Plot Time Series  ------------------------------------

# Plot of the normal TS

plt.plot(AAPL)
plt.title("TS of AAPL")
plt.xlabel("Dates")
plt.ylabel("Values of AAPL")
plt.show()

plt.plot(SPX)
plt.title("TS of SPX")
plt.xlabel("Dates")
plt.ylabel("Values of SPX")
plt.show()

# ------------------------------ Compute Log_returns --------------------------------------

logAAPL = np.diff(np.log(AAPL))
logSPX = np.diff(np.log(SPX))

plt.plot(logAAPL)
plt.plot(logSPX)
plt.title("Log returns")
plt.xlabel("Dates")
plt.ylabel("Log returns")
plt.show()

# ----------------------------------- Regressions -------------------------------------

df1 = sm.add_constant(logSPX)
model = sm.OLS(logAAPL, df1)
resultsRegression = model.fit()

intercept, slope = resultsRegression.params
print(resultsRegression.summary())
print("\n")

# Plotting of the intercept and slope
plt.figure(figsize=(10, 6))
plt.scatter(logSPX, logAAPL, label='Data points')
plt.plot(logSPX, slope * logSPX + intercept, color='red', label='Linear regression')
plt.xlabel('SPX log returns')
plt.ylabel('AAPL log returns')
plt.title('Linear Regression of AAPL log returns on SPX log returns')
plt.legend()
plt.grid(True)
plt.show()

# The slope is the coefficient related to the logSPX (x1) in the model computed
# above, the result is SLOPE = 1.1352. It means that for every change in the logSPX,
# the logAAPL will change 1.1352 times more.

# -------------------------------------- YearFrac -------------------------------------

yf_First_Last = yearfrac(AAPL.index[0], AAPL.index[len(AAPL)-1], 3)

print("The year fraction between the First and Last of AAPL is: ", yf_First_Last, "\n")

# ------------------------------------- Interpolate ------------------------------------

# Import of the data required
xVector = [0, 1, 2, 3, 4]
yVector = [1, 2, 3.5, 4, 2]

xq = 2.7

# Interpolation
yq = np.interp(xq, xVector, yVector)

print("The interpolated valued of ", xq, " is: ", yq, "\n")

# -------------------------------------- Simulation -----------------------------

# Study of the variance
std_variance = []

for i in range(1, 10000, 10):
    rv = np.random.normal(size=i)
    std_variance.append((np.std(rv))**2)

plt.plot(range(1, 10000, 10), std_variance)
plt.xlabel('Number of simulations')
plt.ylabel('Variance of the Gaussian')
plt.title('Asymptotic value of the Std gaussian variance')
plt.grid(True)
plt.show()

# Study of the quantile
rv = np.random.normal(size=100000)
print("The value of the CdF od a std gaussian with q = 0.9 is:", np.quantile(rv, 0.9), "\n")
print("The theoretical value of the 0.9 quantile from the tables is: 1.28155156\n")

# ------------------------------------- Minimization ----------------------------------------------

# The analytical minimization comes from the derivative of first order:
# d/dx = 2*(x-3) = 0 >> x = 3
# d/dy = 2*(y-7) = 0 >> y = 7
# Thus the couple that minimizes the function is (x, y) = (3, 7)

def fxy(xy):
    x, y = xy
    return (x - 3) ** 2 + (y - 7) ** 2


# The numerical minimization comes from the application of the Python functions:
result = sc.optimize.minimize(fxy, [0, 0])

print("The values that minimized the function are X = ", result.x[0], " and Y = ", result.x[1], "\n")
