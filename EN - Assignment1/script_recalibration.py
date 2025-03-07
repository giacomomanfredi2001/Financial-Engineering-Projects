"""
Main script to run the recalibration experiments
"""
# Author: Alessandro Brusaferri
# License: Apache-2.0 license
import numpy as numpy
import scipy.stats as stats
import matplotlib.pyplot as plt

from tools.PrTSF_Recalib_tools import PrTsfRecalibEngine, load_data_model_configs
from tools.prediction_quantiles_tools import plot_quantiles

def calculate_rmse(y_true, y_pred):
    """
    Calculate Root Mean Squared Error (RMSE)
    :param y_true: Array of true values
    :param y_pred: Array of predicted values
    :return: RMSE value
    """
    assert len(y_true) == len(y_pred), "Lengths of true and predicted arrays must be equal."

    squared_diff = (y_true - y_pred) ** 2
    mean_squared_diff = numpy.mean(squared_diff)
    rmse = numpy.sqrt(mean_squared_diff)
    return rmse

def calculate_mae(y_true, y_pred):
    """
    Calculate Mean Absolute Error (MAE)
    :param y_true: Array of true values
    :param y_pred: Array of predicted values
    :return: MAE value
    """
    assert len(y_true) == len(y_pred), "Lengths of true and predicted arrays must be equal."

    abs_diff = numpy.abs(y_true - y_pred)
    mae = numpy.mean(abs_diff)
    return mae

def calculate_smape(y_true, y_pred):
    """
    Calculate Symmetric Mean Absolute Percentage Error (sMAPE)
    :param y_true: Array of true values
    :param y_pred: Array of predicted values
    :return: sMAPE value
    """
    assert len(y_true) == len(y_pred), "Lengths of true and predicted arrays must be equal."

    numerator = numpy.abs(y_true - y_pred)
    denominator = (numpy.abs(y_true) + numpy.abs(y_pred)) / 2

    smape = 100 * numpy.mean(numerator / denominator)
    return smape

def confidence_interval_gaussian(model, y_true, y_pred, alpha=0.05):
    """
    Calculates the confidence interval for RMSE assuming a Gaussian distribution.

    Args:
      y_true (np.ndarray): Ground truth values.
      y_pred (np.ndarray): Predicted values.
      alpha (float, optional): Level of confidence (default: 0.05 for 95% CI).

    Returns:
      tuple: Lower and upper bounds of the confidence interval.
    """
    n = len(y_true)

    error_std = numpy.sqrt(numpy.mean((y_true - y_pred) ** 2) / (n - 2))

    # Determine critical z-value
    z_critical = stats.norm.ppf(1 - alpha / 2)

    # Compute confidence interval for standardized RMSE
    lower_bound = model - z_critical * error_std
    upper_bound = model + z_critical * error_std

    return lower_bound, upper_bound

#--------------------------------------------------------------------------------------------------------------------
# Set PEPF task to execute
PF_task_name = 'EM_price'
# Set Model setup to execute
exper_setup = 'point-ARX'

#---------------------------------------------------------------------------------------------------------------------
# Select run
run_id = 'recalib_opt_grid_1_1'
# Load hyperparams from file (select: load_tuned or optuna_tuner)
hyper_mode = 'load_tuned'
# Plot train history flag
plot_train_history=True
plot_weights=True
#---------------------------------------------------------------------------------------------------------------------

# Load experiments configuration from json file
configs=load_data_model_configs(task_name=PF_task_name, exper_setup=exper_setup, run_id=run_id)

# Instantiate recalibratione engine
PrTSF_eng = PrTsfRecalibEngine(data_configs=configs['data_config'],
                               model_configs=configs['model_config'])

# Get model hyperparameters (previously saved or by tuning)
model_hyperparams = PrTSF_eng.get_model_hyperparams(method=hyper_mode, optuna_m=configs['model_config']['optuna_m'])

# Exec recalib loop over the test_set samples, using the tuned hyperparams
test_predictions = PrTSF_eng.run_recalibration(model_hyperparams=model_hyperparams,
                                               plot_history=plot_train_history,
                                               plot_weights=plot_weights)

# Plot test predictions
plot_quantiles(test_predictions, target=PF_task_name)

print('Done!')

# Evaluation metrics
print('RMSE \n');
rmse_value = calculate_rmse(test_predictions.iloc[:,1], test_predictions.iloc[:,0])
print(rmse_value)
print(confidence_interval_gaussian(rmse_value, test_predictions.iloc[:,1], test_predictions.iloc[:, 0], 0.05))

print('MAE \n');
mae_value = calculate_mae(test_predictions.iloc[:,1], test_predictions.iloc[:,0])
print(mae_value)
print(confidence_interval_gaussian(mae_value, test_predictions.iloc[:,1], test_predictions.iloc[:, 0], 0.05))

print('sMAPE \n');
smape_value = calculate_smape(test_predictions.iloc[:,1], test_predictions.iloc[:,0])
print(smape_value)
print(confidence_interval_gaussian(smape_value, test_predictions.iloc[:,1], test_predictions.iloc[:, 0], 0.05))

#Implementation hourly evalution
predictions = numpy.array(test_predictions.iloc[:, 0])
predictions = predictions.reshape(10, 24)
targets = numpy.array(test_predictions.iloc[:, 1])
targets = targets.reshape(10, 24)

RMSE = numpy.zeros((1, 24))
MAE = numpy.zeros((1, 24))
sMAPE = numpy.zeros((1, 24))

for i in range(24):
    y_true = targets[:, i].transpose()
    y_pred = predictions[:, i].transpose()
    RMSE[0, i] = calculate_rmse(y_true, y_pred)
    MAE[0, i] = calculate_mae(y_true, y_pred)
    sMAPE[0, i] = calculate_smape(y_true, y_pred)

RMSE = RMSE.squeeze()
MAE = MAE.squeeze()
sMAPE = sMAPE.squeeze()

print("Hourly RMSE:", RMSE)
print("Hourly MAE:", MAE)
print("Hourly sMAPE:", sMAPE)

plt.plot(RMSE)
plt.xlabel('Hours')
plt.ylabel('Hourly RMSE')
plt.title('Hourly RMSE')
plt.show()



