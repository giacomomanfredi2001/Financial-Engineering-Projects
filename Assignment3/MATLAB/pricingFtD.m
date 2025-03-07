function spreadCDS = pricingFtD(nSim, rho, dates, discounts, datesCDS, spreadsCDS, recovery)
% Computation of the spread of a First to Default through the Li Model and
% the use of a Gaussian Copula
% 
%INPUT:
% nSim:           number of simulations
% rho:            correlation coefficient
% dates:          vector of dates of the ZR curve
% discounts:      vector of discounts of the ZR curve
% datesCDS:       matrix of reset dates of the CDS (ISP/UCG)
% spreadCDS:      matrix of spread at reset dates of the CDS (ISP/UCG)
% recovery:       vector of the recovery values of ISP/UCG
% 
%USES:
% function zeroRates(dates, discounts)
% function computation_FtD_feeleg(spread, interdelta, discounts, tau_total, indexes_UCG_default, indexes_ISP_default, indexes_no_default, yf_swap)
% function computation_FtD_contingentleg(dates, discounts, tau_total, indexes_UCG_default, indexes_ISP_default, recovery)

    %% Introduction of the quantities

    correlationMatrix = [1 rho; rho 1];
    meanVector = [0; 0];

    % Divisions of the matrixes
    datesCDS_ISP = datesCDS;
    datesCDS_UCG = datesCDS;

    spreadsCDS_ISP = spreadsCDS(:, 1);
    spreadsCDS_UCG = spreadsCDS(:, 2);

    recovery_ISP = recovery(1);
    recovery_UCG = recovery(2);

    % Standard creation
    std_act_365 = 3;
    std_30_360 = 6;
    
    %% 1st step: Computation of the ui
    
    % Computation of the multivariate samples
    rng(0);
    samplesX = mvnrnd(meanVector,correlationMatrix,nSim);
    
    u = normcdf(samplesX);
    u_UCG = u(:, 1);
    u_ISP = u(:, 2);
    
    %% 2nd step: Inverting the structure of ui to get time to default
    
    [~, survProbs_exact_ISP, intensities_ISP] = bootstrapCDS(dates, discounts, datesCDS_ISP, spreadsCDS_ISP, 2, recovery_ISP);
    [~, survProbs_exact_UCG, intensities_UCG] = bootstrapCDS(dates, discounts, datesCDS_UCG, spreadsCDS_UCG, 2, recovery_UCG);
    
    % Modification of the datesCDS
    yfCDS = [0; yearfrac(dates(1), datesCDS, std_act_365)];
    
    % We insert a prob = 0 at the end of the survProb to always obtain an index
    survProbs_exact_ISP = [survProbs_exact_ISP; 0];
    survProbs_exact_UCG = [survProbs_exact_UCG; 0];
    
    % UCG and ISP inverting
    tau_UCG = zeros(nSim, 1);
    tau_ISP = zeros(nSim, 1);
    
    for i = 1: nSim
    
        % UCG
        index_UCG = find(survProbs_exact_UCG <= u_UCG(i), 1);
    
        if(index_UCG <= 5)
            tau_UCG(i) = - (1/intensities_UCG(index_UCG - 1)) * log(u_UCG(i)/survProbs_exact_UCG(index_UCG-1)) + yfCDS(index_UCG-1);
        else
            tau_UCG(i) = 10;      % Placeholder value
        end
        
    
        % ISP
        index_ISP = find(survProbs_exact_ISP <= u_ISP(i), 1);
    
        if(index_ISP <= 5)
            tau_ISP(i) = - (1/intensities_ISP(index_ISP - 1)) * log(u_ISP(i)/survProbs_exact_ISP(index_ISP-1)) + yfCDS(index_ISP-1);
        else
            tau_ISP(i) = 10;      % Placeholder value
        end
    end
    
    %% 3rd step: Find the defaultable obligor
    tau_total = [tau_ISP tau_UCG];
    
    % Find those defaultable ISP, UCG, no default
    indexes_UCG_default = find(tau_UCG < tau_ISP);
    indexes_ISP_default = find(tau_UCG > tau_ISP);
    indexes_no_default = find(tau_UCG == tau_ISP);
    
    %% 4th step: Spread computation
    
    % Interdelta
    delta_30_360 = [0; yearfrac(dates(1), datesCDS(1:4), std_30_360)];
    interdelta = delta_30_360(2:end) - delta_30_360(1:end-1);
    
    % Discounts
    delta_ACT_365_total = yearfrac(dates(1), dates(2:end), 3);
    zero_rates = zeroRates(dates, discounts)/100;
    zr_1y_swap = interp1(delta_ACT_365_total, zero_rates(2:end), 1, 'linear');
    discount_1y_swap = exp(-zr_1y_swap*1);
    
    discounts_fee = [discount_1y_swap; discounts(12:14)];
    
    yfSwap = yearfrac(dates(1), datesCDS(1:4), std_30_360);

    % Computation of the spread
    options = optimset('Display', 'Off');
    spreadCDS = lsqnonlin(@(s) computation_FtD_feeleg(s, interdelta, discounts_fee, tau_total, indexes_UCG_default, indexes_ISP_default, indexes_no_default, yfSwap) - ...
        computation_FtD_contingentleg(dates, discounts, tau_total, indexes_UCG_default, indexes_ISP_default, [recovery_ISP; recovery_UCG]), 0, 0, Inf, options);

end % function pricingFtD