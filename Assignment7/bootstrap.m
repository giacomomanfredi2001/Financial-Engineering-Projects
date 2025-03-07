function [dates, discounts]=bootstrap(datesSet, ratesSet)
% Bootstrap function to calculate the discount values curve
% 
%INPUT
% datesSet:     structure containing dates of Depo, Futures and Swaps
% ratesSet:     structure containing rates of Depo, Futures and Swaps
% 
%OUTPUT
% dates:        values of dates interpolated
% discounts:    values of rates interpolated
% 
%USES
% function interpolation(B1, B2, yf1, yf2, yf)

%% Quantities of interests

% Division of dates
settlementDate = datesSet.settlement;            % Date of settlement

deposDates = datesSet.depos;                     % Dates of IB deposits
futuresDates = datesSet.futures(1:7, :);         % Dates of Futures
swapsDates = datesSet.swaps(1:end);              % Dates of Swaps

% Division of rates - LIBOR

deposRates = ratesSet.depos(:, :);               % Rates of IB deposits
futuresRates = ratesSet.futures(1:7, :);         % Rates of Futures
swapsRates = ratesSet.swaps(1:end, :);           % Rates of Swaps

%% IB deposit
% The standard for the IB deposit is the ACT/360

% Compute the index of the first depo after the first future
index = find(deposDates >= futuresDates(1, 1), 1);
deposDates = deposDates(1:index);
deposRates = deposRates(1:index, :);

% Compute the mid LIBOR for IB and year fraction
mid_L_IB = mean(deposRates, 2);
yf_IB = yearfrac(settlementDate, deposDates, 2);  

% Initialize Discount fixed vectors
B_final = zeros(1, index + 14 + 1);           
B_final_to_give_back = zeros(1, index +7 + 1);

B_final(1:index) = 1./(1 + yf_IB.*mid_L_IB); 
B_final_to_give_back(1:index) = 1./(1 + yf_IB.*mid_L_IB); 

% Up to now, we've computed B(t0, ti) only for the IB deposit

%% Futures

% We need to consider only the first 7 Futures with  ACT/360

% Computation of the Mid Prices and year fraction of the futures
mid_L_F = mean(futuresRates, 2);
yf_IB = yearfrac(settlementDate, deposDates, 3);
yf_F = yearfrac(settlementDate, futuresDates(:, 1), 3);  

yf_IB_F = yf_IB;                   % delta(t0, ti) [DEPO FUTURE]

% Helping variable to fill the vector, reminding at which step I'm up to
% now
count = 1;

for i = index+1:2:index+14

    % Add the settlement date year fraction
    yf_IB_F = [yf_IB_F; yf_F(count)];    

    if (yf_IB_F(i) <= yf_IB_F(i-1))            

        % INTERPOLATION
        B_temp = interpolation(B_final(i-2), B_final(i-1), ...
            yf_IB_F(i-2), yf_IB_F(i-1), yf_IB_F(i));

        B_supp = B_final(i-1);
        B_final(i-1) = B_temp;
        B_final(i) = B_supp;

        % Invert yf to the correct order
        D_temp = yf_IB_F(i);
        yf_IB_F(i) = yf_IB_F(i-1);
        yf_IB_F(i-1) = D_temp;
        

    else                                            
        
        % EXTRAPOLATION FLAT
        % Introduction of a flat interpolation on the rates due to the 
        % guide lines given by the professor

        zeroRate_temp = - log(B_final(i-1))/yf_IB_F(end-1);
        
        B_temp = exp(-zeroRate_temp*yf_IB_F(end));
        B_final(i) = B_temp;

    end

    % Now compute the 2nd part: B(t0; ti, ti+1)
    dti = yearfrac(futuresDates(count, 1), futuresDates(count, 2), 2);

    B_temp = B_temp / (1 + dti*mid_L_F(count));

    % Update the final values for the computation

    B_final(i+1) = B_temp;                            % UPDATED B(t0, ti) for settlement and expiry vector
    B_final_to_give_back(index + count) = B_temp;     % UPDATED B(t0, ti) for expiry vector

    yf_IB_F = [ yf_IB_F; yearfrac(settlementDate, futuresDates(count, 2), 3)];

    count = count +1;

end

%% Swaps
% The standard to consider is the 30/360 (European)

yf_S = yearfrac(settlementDate, swapsDates(1), 3);

yf_IB_F_S = [yf_IB_F; yf_S];

% INTERPOLATION

index = find(yf_IB_F_S >= yf_IB_F_S(end), 1);
B_temp = interpolation(B_final(index-1), B_final(index), ...
        yf_IB_F_S(index-1), yf_IB_F_S(index), yf_IB_F_S(end));

% Update the final value
B_final_to_give_back(end) = B_temp;               % UPDATED B(t0, ti) only expiry
B_final(end) = B_temp;                            % UPDATED B(t0, ti) settement/expiry

%% Swaps - COMPLETENESS
% The standard to consider is the 30/360 (European)

% Computation of mid rates
mid_L_S = mean(swapsRates, 2);
swapsDates_complete = [settlementDate; swapsDates];

yf_complete_S_fixed = yearfrac(swapsDates_complete(1:end-1), swapsDates_complete(2:end), 6);

% Initialization of the Discount vector of the Swap
B_swap = zeros(size(swapsRates, 1), 1);
B_swap(1) = B_final(end);

% Computation of the Complete Swap curve
for i=2:length(B_swap) 
    B_swap(i) = (1 - mid_L_S(i)*(yf_complete_S_fixed'*B_swap))...
        /(1 + yf_complete_S_fixed(i)*mid_L_S(i));
end

%% Creation of the return vectors

% Discount factors
discounts = [1; B_final_to_give_back(1:end-1)' ; B_swap(2:end)];

% Dates
dates = [settlementDate; deposDates ; futuresDates(:, 2); swapsDates(2:end)];


end %function bootstrap