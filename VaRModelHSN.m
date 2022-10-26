clear all; close all; clc;
load DataSTDoil.mat
%% In this section, we'll go through various techniques for calculating VAR and relative backtesting using the following steps:
%  1. Normal Distribution approach (Parametric)
%  2. Historical Simulation approach (Non-Parametric)
%  3. Extreme Value Theory approach (Semi-Parametric)

%% Arrange data to make the analysis that follows more fluid.
% Define the window
oil_date = oil_date.Date;
TestWindowStart = find(year(oil_date)==2019,1);
TestWindowEnd = find(year(oil_date)==2022,1,'last');
TestWindow = TestWindowStart:TestWindowEnd;
WindowSize = 250;
index = oil_std_returns.Properties.VariableNames;
% Define the VaR confindance level
pVaR = [0.05,0.01];

%% Normal distribution approach using movable windows (Paramatric)
Zscore = norminv(pVaR);

% Calculating the VaR in a for loop for every point of test window and for
% every companies
% NB: the var is calculated for the negative side of the distribution of
% residuals
Normal95 = zeros(length(TestWindow),width(index));
Normal99 = zeros(length(TestWindow),width(index));

for j = 1:width(index)
    for t=TestWindow
        i=t-TestWindowStart+1;
        EstimationWindow = t-WindowSize:t-1;
        Sigma_window = std(oil_std_returns{EstimationWindow,index(j)});
        Normal95(i,j) = -Zscore(1)*Sigma_window;
        Normal99(i,j) = -Zscore(2)*Sigma_window;
    end

% Rappresent the VaR for each confidance on a figure
    figure('Name','VaR Normal Distribution Method')
    plot(oil_date(TestWindow),[Normal95(:,j) Normal99(:,j)])
    xlabel('Date');
    ylabel('VaR')
    legend({'95% Confidence level','99% Confidence level'},'location','Best');
    title(strcat(['VaR Normal Distribution method '],index(j)));

% VaR Backtesting of Normal distribution approach
    Returns_test = oil_std_returns{TestWindow,index(j)};
    Date_test = oil_date(TestWindow);
    figure('Name','Return Test and VaR Estimation')
    plot(Date_test,[Returns_test -Normal95(:,j) -Normal99(:,j)]);
    xlabel('VaR Estimated');
    ylabel('Date')
    title(strcat(['Returns test vs VaR '],index(j)));
    legend({'Return Test','Normal VaR 95','Normal VaR 99'});
    vbt = varbacktest(Returns_test,[Normal95(:,j) Normal99(:,j)],'PortfolioID',index(j),'VaRID',{'Normal95','Normal99'},'VaRLevel',0.95);
    backtest = runtests(vbt);
    disp(backtest)
end

%% Historical Simulation Approach using movable windows
    Historical95 = zeros(length(TestWindow),width(index));
    Historical99 = zeros(length(TestWindow),width(index));

for j = 1:width(oil_std_returns)
    for t = TestWindow
        i = t-TestWindowStart+1;
        EstimationWindow = t-WindowSize:t-1;
        Historical95(i,j) = -quantile(oil_std_returns{EstimationWindow,index(j)},pVaR(1));
        Historical99(i,j) = -quantile(oil_std_returns{EstimationWindow,index(j)},pVaR(2));
    end
    figure('Name','VaR Historical Simulation Method')
    plot(oil_date(TestWindow),[Historical95(:,j),Historical99(:,j)]);
    xlabel('VaR Historical')
    ylabel('Date')
    title(strcat(['VaR Historical Method'],index(j)));
    legend({'Historical 95','Historical 99'});

    % VaR backtesting for Historical simulation approach
    Returns_test = oil_std_returns{TestWindow,index(j)};
    Date_test = oil_date(TestWindow);
    figure('Name','Return Test and Historical VaR')
    plot(Date_test,[Returns_test,-Historical95(:,j),-Historical99(:,j)]);
    xlabel('Date')
    ylabel('Returns Test and Historical VaR');
    title(strcat(['VaR and Returns Test for'],index(j)));
    legend({'Returns','Historical 95','Historical 99'},'location','best');
    vbt =varbacktest(Returns_test,[Historical95(:,j),Historical99(:,j)],'PortfolioID',index(j),'VaRID',{'Historical 95','Historical 99'},'VaRLevel',1-pVaR(1));
    backtest = runtests(vbt);
    disp(backtest)
end

%% Building graphs in order to compare the two method at 95 % of confidance on the whole period
for j = 1:width(index)
    VaRData = [-Normal95(:,j),-Historical95(:,j)];
    VaRFormat = {'-','--'};
    IndexNormal95 = oil_std_returns{TestWindow,index(j)}<=VaRData(:,1); %% Da vedere
    IndexHistorical95 = oil_std_returns{TestWindow,index(j)}<=VaRData(:,2);
    figure('Name','VaR Excedance')
    bar(oil_date(TestWindow),oil_std_returns{TestWindow,index(j)},'FaceColor',['k',0.7,0.7]);
    hold on
    for i=1:width(VaRData)
        stairs(oil_date(TestWindow),VaRData(:,i),VaRFormat{i});
    end
    xlabel('Date')
    ylabel('VaR Normal and Historical with Exceedances')
    title(strcat(['VaR 95% violations for Normal and Historical of'],index(j)));
    ax = gca;
    ax.ColorOrderIndex = 2;
    plot(oil_date(TestWindow(IndexNormal95)),-Normal95(IndexNormal95,j),'o',oil_date(TestWindow(IndexHistorical95)),-Historical95(IndexHistorical95,j),'x','MarkerSize',8,'LineWidth',1);
    hold off
    legend({'Returns','Normal','Historical','Normal Violations','Historical Violations'},'Location','best')
end


%% Conclusions for Normal Distribution
   % PortfolioID      VaRID       VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    __________    ________    _____    ______    ______    ______    ______    ______    ______    ______

     %"CVX_res"     "Normal95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
     %"CVX_res"     "Normal99"      0.95      green    reject    reject    accept    reject    accept    reject    reject

   % PortfolioID      VaRID       VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    __________    ________    _____    ______    ______    ______    ______    ______    ______    ______

    % "COP_res"     "Normal95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
    % "COP_res"     "Normal99"      0.95      green    reject    reject    accept    reject    accept    reject    reject

   % PortfolioID      VaRID       VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    __________    ________    _____    ______    ______    ______    ______    ______    ______    ______

    % "TOT_res"     "Normal95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
    % "TOT_res"     "Normal99"      0.95      green    reject    reject    accept    reject    accept    reject    reject

   % PortfolioID      VaRID       VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    __________    ________    _____    ______    ______    ______    ______    ______    ______    ______

    % "ENI_res"     "Normal95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
    % "ENI_res"     "Normal99"      0.95      green    reject    reject    accept    reject    accept    reject    reject
    
%% Conclusions for Historical Simulation
    %PortfolioID         VaRID         VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
    %___________    _______________    ________    _____    ______    ______    ______    ______    ______    ______    ______

     %"CVX_res"     "Historical 95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
     %"CVX_res"     "Historical 99"      0.95      green    reject    reject    accept    reject    accept    reject    reject

    %PortfolioID         VaRID         VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
    %___________    _______________    ________    _____    ______    ______    ______    ______    ______    ______    ______

     %"COP_res"     "Historical 95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
     %"COP_res"     "Historical 99"      0.95      green    reject    reject    reject    reject    accept    reject    reject

   % PortfolioID         VaRID         VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    _______________    ________    _____    ______    ______    ______    ______    ______    ______    ______

    % "TOT_res"     "Historical 95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
     %"TOT_res"     "Historical 99"      0.95      green    reject    reject    accept    reject    accept    reject    reject

   % PortfolioID         VaRID         VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
    %___________    _______________    ________    _____    ______    ______    ______    ______    ______    ______    ______

     %"ENI_res"     "Historical 95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
     %"ENI_res"     "Historical 99"      0.95      green    reject    reject    reject    reject    accept    reject    reject















