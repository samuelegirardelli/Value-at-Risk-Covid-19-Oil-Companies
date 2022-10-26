%% calculating the financial series' returns for the oil firm, in particularly CVX COP TOT ENI, from 1-Jan-2018 to 30-Aug-2022

clear all; close all; clc;

% Loading the database from Excel. The closing price has been downloaded from
% Datasource.
cvx_price = table2timetable(readtable('Database.xlsx','Sheet','CVX'));
cvx_price = renamevars(cvx_price,'Close','CVX');
cop_price = table2timetable(readtable('Database.xlsx','Sheet','COP'));
cop_price = renamevars(cop_price,'Close','COP');
tot_price = table2timetable(readtable('Database.xlsx','Sheet','TOT'));
tot_price = renamevars(tot_price,'Close','TOT');
eni_price = table2timetable(readtable('Database.xlsx','Sheet','ENI'));
eni_price = renamevars(eni_price,'Close','ENI');

oil_prices = synchronize(cvx_price,cop_price,tot_price,eni_price,'first');
oil_prices = fillmissing(oil_prices,'nearest');


% Calculating the log returns
cvx_return = timetable(cvx_price.Date(2:end),tick2ret(cvx_price.CVX,'Method','continuous'),'VariableNames',{'CVX'});
cop_return = timetable(cop_price.Date(2:end),tick2ret(cop_price.COP,'method','continuous'),'VariableNames',{'COP'});
tot_return = timetable(tot_price.Date(2:end),tick2ret(tot_price.TOT,'method','continuous'),'VariableNames',{'TOT'});
eni_return = timetable(eni_price.Date(2:end),tick2ret(eni_price.ENI,'method','continuous'),'VariableNames',{'ENI'});

oil_returns = synchronize(cvx_return,cop_return,tot_return,eni_return,'first');
oil_returns = fillmissing(oil_returns,'nearest');


% Clearing the workspace from useful data
clear('cop_price','cop_return','cvx_price','cvx_return','eni_price','eni_return','tot_price','tot_return');

% Saving the file
save OilReturns

%% Visualization of the Oil Company's data

% Closing Price plot with base 100 using plot tool
figure('Name','Closing prices of oil companies')
hold on
grid on
plot(oil_prices.Date,ret2price(oil_returns.CVX)*100,'Color','red');
plot(oil_prices.Date,ret2price(oil_returns.COP)*100,'Color','green');
plot(oil_prices.Date,ret2price(oil_returns.TOT)*100,'Color','cyan');
plot(oil_prices.Date,ret2price(oil_returns.ENI)*100,'Color','yellow');

% Font
datetick('x')
xlabel('Date','FontSize',12,'FontWeight','bold');
ylabel('Companies Price Values','FontSize',12,'FontWeight','bold');
title('Relative Daily Closing Prices')
legend('CVX','COP','TOT','ENI','location','best');

hold off

%% Graphing the results
figure('Name','Returns of oil companies')
%CVX Return
subplot(2,2,1);
plot(oil_returns.Time,oil_returns.CVX)
xlabel('Date','FontSize',8,'FontWeight','bold');
ylabel('CVX Return','FontSize',8,'FontWeight','bold');
%COP Return
subplot(2,2,2);
plot(oil_returns.Time,oil_returns.COP)
xlabel('Date','FontSize',8,'FontWeight','bold');
ylabel('COP Return','FontSize',8,'FontWeight','bold');
%TOT Return
subplot(2,2,3);
plot(oil_returns.Time,oil_returns.TOT)
xlabel('Date','FontSize',8,'FontWeight','bold');
ylabel('TOT Return','FontSize',8,'FontWeight','bold');
%ENI Return
subplot(2,2,4);
plot(oil_returns.Time,oil_returns.ENI)
xlabel('Date','FontSize',8,'FontWeight','bold');
ylabel('ENI Return','FontSize',8,'FontWeight','bold');

%% Looking if some empirical distributions fit our data - Normal
figure('Name','t-Student distribution fitting')
subplot(2,2,1)
histfit(oil_returns.CVX,'Normal');
subplot(2,2,2)
histfit(oil_returns.COP,'Normal');
subplot(2,2,3)
histfit(oil_returns.TOT,'Normal');
subplot(2,2,4)
histfit(oil_returns.ENI,'Normal');

%% Validanting the assumption of independence and identically distributed

%% Autocorrelation function verification
figure('Name','ACF graphs of oil companies');
% CVX ACF
subplot(2,2,1)
autocorr(oil_returns.CVX);
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('CVX ACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])
% COP ACF
subplot(2,2,2)
autocorr(oil_returns.COP);
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('COP ACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])
% TOT ACF
subplot(2,2,3)
autocorr(oil_returns.TOT);
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('TOT ACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])
% ENI ACF
subplot(2,2,4)
autocorr(oil_returns.ENI);
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('ENI ACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])

%% Partial autocorrelation function verification

figure('Name','PACF for oil companies')
% CVX PACF
subplot(2,2,1)
parcorr(oil_returns.CVX)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('CVX PACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])
% COP PACF
subplot(2,2,2)
parcorr(oil_returns.COP)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('COP PACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])
% TOT PACF
subplot(2,2,3)
parcorr(oil_returns.TOT)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('TOT PACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])
% ENI PACF
subplot(2,2,4)
parcorr(oil_returns.ENI)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('ENI PACF','FontSize',8,'FontWeight','bold')
ylim([-0.2 1])

%% Checking the dependency of the square residuals
figure('Name','Squared Residuals ACF')
% CVX ACF Squared Residuals
subplot(2,2,1)
autocorr(oil_returns.CVX.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('CVX ACF Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])
% COP ACF Squared Residuals
subplot(2,2,2)
autocorr(oil_returns.COP.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('COP ACF Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])
% TOT ACF Squared Residuals
subplot(2,2,3)
autocorr(oil_returns.TOT.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('TOT ACF Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])
% ENI ACF Squared Residuals
subplot(2,2,4)
autocorr(oil_returns.ENI.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('ENI ACF Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])

%% Working with IID data: we must normalize our returns in order to do accurate analysis.
% Summary of ours data
% CVX: 1217×1 double     COP: 1217×1 double     TOT: 1217×1 double     ENI: 1217×1 double
% Values:              % Values:              % Values:              % Values:
% Min         -0.25006 % Min       -0.28555   % Min       -0.24116   % Min         -0.23385
% Median    0.00025825 % Median           0   % Median           0   % Median    0.00029176
% Max           0.2049 % Max        0.22485   % Max         0.1308   % Max          0.13916

% Modelling the first moment of the distribution MEAN and the
% second moment VARIANCE for the whole dataset.

% While fitting a model for the mean, we notice that our data appears to
% already be devoid of a mean from a first study of the summary.

% ARIMA Model and GARCH Model togheter 
model = arima('AR',NaN,'Distribution','Gaussian','Variance',gjr(1,1));
option = optimoptions(@fmincon,'Display','off','Diagnostic','off','Algorithm','sqp','TolCon',1e-7);

% CVX Estimate and inference
cvx_fit = estimate(model,oil_returns.CVX,'Option',option);
[~,CVX_variances] = infer(cvx_fit,oil_returns.CVX);

% COP Estimate and inference
cop_fit = estimate(model,oil_returns.COP,'Option',option);
[~,COP_variances] = infer(cop_fit,oil_returns.COP);

% TOT Estimate and inference
tot_fit = estimate(model,oil_returns.TOT,'Option',option);
[~,TOT_variances] = infer(tot_fit,oil_returns.TOT);

% ENI Estimate and inference
eni_fit = estimate(model,oil_returns.ENI,'Option',option);
[~,ENI_variances] = infer(eni_fit,oil_returns.ENI);
%% Remarks
% We don't derive any significant conclusions about the AR process from the 
% model that we create.
% Conversly the whole sample seems to have significative variance
% dependence among time and the parameter of the GARCH model result
% significant for at least 5%.
% Some significance has been found for the 'leverage' GARCH parameter.
%% Creating dataframe with standardized returns in order to apply analysis.

% CVX Standardized Residuals and look for confirmation
cvx_std = oil_returns.CVX./sqrt(CVX_variances);

figure('Name','ACF of Squared STD Residuals')

subplot(2,2,1)
autocorr(cvx_std.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('CVX ACF STD Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])

% COP Standardized Residuals and look for confirmation
cop_std = oil_returns.COP./sqrt(COP_variances);

subplot(2,2,2)
autocorr(cop_std.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('COP ACF STD Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])

% TOT Standardized Residuals and look for confirmation
tot_std = oil_returns.TOT./sqrt(TOT_variances);

subplot(2,2,3)
autocorr(tot_std.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('TOT ACF STD Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])

% ENI Standardized Residuals and look for confirmation
eni_std = oil_returns.ENI./sqrt(ENI_variances);

subplot(2,2,4)
autocorr(eni_std.^2)
xlabel('Lag','FontSize',8,'FontWeight','bold')
ylabel('ENI ACF STD Squared Residuals','FontSize',8,'FontWeight','bold')
ylim([0 1])

%% Results: now that we are working on approximately I.I.D. Data, it is confirmed.

% Rearrange the data for the next step of the analysis
oil_std_returns = table(cvx_std,cop_std,tot_std,eni_std,'VariableNames',{'CVX_res','COP_res','TOT_res','ENI_res'});
oil_variances = table(CVX_variances,COP_variances,TOT_variances,ENI_variances,'VariableNames',{'CVX_var','COP_var','TOT_var','ENI_var'});
oil_ret = table(oil_returns.CVX,oil_returns.COP,oil_returns.TOT,oil_returns.ENI,'VariableNames',{'CVX_ret','COP_ret','TOT_ret','ENI_ret'});
oil_date = table(oil_returns.Time,'VariableNames',{'Date'});
clear('COP_residuals','cop_std','COP_variances','CVX_residuals','cvx_std','CVX_variances','ENI_residuals','eni_std','ENI_variances','model','option','TOT_residuals','tot_std','TOT_variances','oil_prices','oil_returns');
save DataSTDoil.mat


















