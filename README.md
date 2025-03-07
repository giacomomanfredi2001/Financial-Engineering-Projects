## Financial-Engineering-Projects

This repository is a collection of the projects and assignments of the course of Financial Engineering (2024).
The projects cover a broad range of topics related to financial derivatives, pricing models, numerical methods, risk management, and machine learning applications. Below is a summary of the key topics covered in each assignment.

# Assignment 1: 
Option Pricing and Error Analysis
European Call Option Pricing: Implemented Black’s formula, the CRR binomial tree, and Monte Carlo (MC) simulations.
Error Rescaling: Evaluated pricing errors across different methods and optimized parameters for improved accuracy.
Exotic Options: Priced a Knock-In Call option using closed-form solutions, CRR, and MC methods.
Vega Sensitivity: Compared numerical and analytical calculations of Vega for various options.
Bermudan Option Pricing: Analyzed the pricing of Bermudan options using the CRR method and compared their behavior with European options.
Monte Carlo Enhancements: Applied variance reduction techniques (antithetic variables) in MC simulations.

# Assignment 2: Yield Curves and Sensitivities
Yield Curve Construction: Developed a discount factor curve and zero rates using deposits, futures, and swaps.
Sensitivity Analysis: Computed DV01, BPV, and duration for interest rate swaps to measure portfolio exposure to rate shifts.
Theoretical Models: Derived bond pricing formulas and examined the Garman-Kohlhagen model for European call options with time-dependent rates and volatility.

# Assignment 3: Asset Swaps, CDS, and Time-Series Analysis
Asset Swap Spread: Calculated spreads over Euribor 3m using market data.
CDS Curve Bootstrapping: Constructed CDS term structures for multiple obligors (e.g., Intesa Sanpaolo) using spline interpolation.
First-to-Default Pricing: Priced First-to-Default (FtD) contracts via the Li model and Monte Carlo simulations.
Python Time-Series Analysis: Implemented basic time-series analytics, including log-return plotting and regression analysis.

# Assignment 4: Value at Risk (VaR) and Expected Shortfall (ES)
Variance-Covariance Method: Estimated VaR and ES at a 99% confidence level using a t-distribution.
Historical Simulation & Bootstrap: Applied historical simulation techniques for risk estimation.
Principal Component Analysis (PCA): Used PCA for dimensionality reduction in VaR calculations.
Monte Carlo VaR: Developed full Monte Carlo and Delta-Normal VaR methodologies.
Cliquet Option Pricing: Investigated the impact of counterparty risk on Cliquet options.

# Assignment 5: Advanced Derivative Pricing and Volatility Surface Calibration
Certificate Pricing: Priced structured products based on ENEL and AXA stocks using Monte Carlo simulations.
Digital Option Pricing: Compared pricing methods (Black-Scholes vs. implied volatility).
Fourier Transform Methods: Applied Lewis' formula and FFT for efficient option pricing.
Volatility Surface Calibration: Used a mean-variance mixture model to calibrate volatility surfaces and analyzed skewness and Vega sensitivity.

# Assignment 6: Interest Rate Risk and Hedging
Yield Curve Extension: Extended bootstrapping techniques to model discount factors beyond 12-year maturities.
Caplet Pricing & Volatility Calibration: Estimated spot volatilities from market cap prices using the Bachelier model.
Structured Product Pricing: Calculated upfront payments for structured interest rate instruments.
Delta & Vega Hedging: Designed hedging strategies using swaps and caps.

# Assignment 7: Bermudan Swaptions and Certificate Pricing
Bermudan Swaption Pricing: Implemented pricing using the Hull-White model via tree-based methods and closed-form solutions.
Certificate Pricing with NIG Model: Priced structured bonds using the Normal Inverse Gaussian (NIG) model, FFT, Variance Gamma, and Monte Carlo simulations.
Black Model Adjustments: Assessed the impact of digital risk adjustments on certificate pricing.
Energy Price and Load Forecasting (EPLF Assignments)

# EPLF 1: Regularization in Energy Price Forecasting
Lasso, Ridge, and Elastic Net Regression: Applied regularization techniques to enhance electricity price forecasts.
Feature Selection: Used LASSO for optimal variable selection in time-series forecasting.
Seasonality Analysis: Examined seasonal effects on energy prices using the Elastic Net model.

# EPLF 2: Deep Neural Networks (DNN) Hyperparameter Tuning
Hyperparameter Optimization: Used Optuna for learning rate and hidden layer size tuning.
Loss Function Analysis: Compared configurations to optimize loss minimization.
Model Performance: Evaluated the effects of overfitting and generalization in DNNs.

# EPLF 3: Distributional Neural Networks & Quantile Regression
Quantile Regression Networks: Used pinball loss functions to estimate data distribution boundaries.
Distributional Neural Networks (DNN): Implemented Normal and Johnson’s SU distributions for probabilistic forecasting.
Model Performance Comparison: Assessed models via Pinball and Winkler scores, with JSU-DNN performing best.
Risk Management (RM Assignments)

# RM 1: Hazard Rate and Z-Spread Calculation
Hazard Rate Curve Bootstrapping: Estimated hazard rates for Investment Grade (IG) and High Yield (HY) bonds.
Z-Spread Calculation: Determined Z-spreads by shifting the zero-rate curve to match defaultable and risk-free bond prices.
Market-Implied Transition Matrices: Constructed rating migration matrices from hazard rates.

# RM 2: Present Value and Credit VaR
Monte Carlo PV Calculation: Simulated present value considering default probabilities.
Credit VaR Simulation: Modeled default correlation effects in Credit VaR calculations.
Concentration Risk Analysis: Evaluated the impact of portfolio diversification on risk exposure.
