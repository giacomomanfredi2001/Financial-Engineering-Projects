function [yf_ACT360, yf_ACT365, yf_30360] = contracts_yf(settlement_date, dates)
% The function allow the multiple computation through the year frac MATLAB
% function writing the 3 most used conventions
% 
% INPUT:
% settlement_date:             date of settlement of the yf to compute
% dates:                       list of dates to compute the yf on
% 
% OUTPUT:
% yf_ACT360:                   yf_computed in ACT/360
% yf_3ACT65:                   yf_computed in ACT/365
% yf_30360:                    yf_computed in 30/360 (EU)

    %% Conventions

    conv_ACT360 = 2;
    conv_ACT365 = 3;
    conv_30360 = 6;

    yf_ACT360 = yearfrac(settlement_date, dates, conv_ACT360);
    yf_ACT365 = yearfrac(settlement_date, dates, conv_ACT365);
    yf_30360 = yearfrac(settlement_date, dates, conv_30360);

end % function contracts_yf