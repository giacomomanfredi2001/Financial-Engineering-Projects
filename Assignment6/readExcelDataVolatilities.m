function [flat_vols,strikes] = readExcelDataVolatilities(filename)
% Reads data from excel
%  It reads flat volatilities and strikes from the quotes
%  All input rates are in % units
%
% INPUTS:
%  filename: excel file name where data are stored
% 
% OUTPUTS:
%  flat_vols: flat volatilities
%  strikes: strikes in the quotes

flat_vols = xlsread(filename, 1, 'F2:R17');
flat_vols = flat_vols * 1e-4;
flat_vols = [flat_vols(1,:);flat_vols(3:end,:)];

strikes = xlsread("Caps_vol_20-2-24.xlsx",1,'F1:R1');
strikes = strikes'/100;

end % function readExcelDataVolatilities