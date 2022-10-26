clear all; close all; clc;
load DataSTDoil.mat

%% We want to analyze our companies with a semi-parametric method: Extreme Value Theory
% Does the EVT perform better than the previous method used to calculate
% the VaR?
% We assume that a 10% tail cutoff will be used.

%% First of all we want to compute an EVT distribution to approximate
n_point = 200;
tail_fr = 0.10;


%% Working on the whole dataset
tails_cvx = paretotails(oil_std_returns.CVX_res,tail_fr,1-tail_fr,'kernel');
tails_cop = paretotails(oil_std_returns.COP_res,tail_fr,1-tail_fr,'kernel');
tails_tot = paretotails(oil_std_returns.TOT_res,tail_fr,1-tail_fr,'kernel');
tails_eni = paretotails(oil_std_returns.ENI_res,tail_fr,1-tail_fr,'kernel');

%% Plotting the distribution using the kernel distribution and the power law decay on the tails of our distributions
figure('Name','GPD on Tails Distributions')
% CVX
subplot(2,2,1)
hold on
grid on
min_prob = cdf(tails_cvx,min(oil_std_returns.CVX_res));
max_prob = cdf(tails_cvx,max(oil_std_returns.CVX_res));

plow_tail_cvx = linspace(min_prob,tail_fr,n_point);
pup_tail_cvx = linspace(max_prob,1-tail_fr,n_point);
pint_cvx = linspace(tail_fr,1-tail_fr,n_point);

plot(icdf(tails_cvx,plow_tail_cvx),plow_tail_cvx,'red','LineWidth',4);
plot(icdf(tails_cvx,pint_cvx),pint_cvx,'black','LineWidth',2);
plot(icdf(tails_cvx,pup_tail_cvx),pup_tail_cvx,'blue','LineWidth',2);

xlabel('Standardized Residuals')
ylabel('Probability')
title('Empirical CDF CVX')
legend({'Pareto Lower Tail','Kernel Smoothed Interior','Pareto Upper Tail'},'location','northwest')
hold off

% COP 
subplot(2,2,2)
hold on 
grid on

min_prob = cdf(tails_cop,min(oil_std_returns.COP_res));
max_prob = cdf(tails_cop,max(oil_std_returns.COP_res));

plow_tail_cop = linspace(min_prob,tail_fr,n_point);
pint_cop = linspace(tail_fr,1-tail_fr,n_point);
pup_tail_cop = linspace(1-tail_fr,max_prob,n_point);

plot(icdf(tails_cop,plow_tail_cop),plow_tail_cop,'red','LineWidth',4);
plot(icdf(tails_cop,pint_cop),pint_cop,'black','LineWidth',2);
plot(icdf(tails_cop,pup_tail_cop),pup_tail_cop,'blue','LineWidth',2);

xlabel('Standardized Residuals')
ylabel('Probability')
title('Empirical Distribution Function COP')
legend({'Pareto Lower Tail','Kernel Smoothed Interior','Pareto Upper Tail'},'location','northwest')
hold off

% TOT
subplot(2,2,3)
hold on 
grid on

min_prob = cdf(tails_tot,min(oil_std_returns.TOT_res));
max_prob = cdf(tails_tot,max(oil_std_returns.TOT_res));

plow_tail_tot = linspace(min_prob,tail_fr,n_point);
pint_tot = linspace(tail_fr,1-tail_fr,n_point);
pup_tail_tot = linspace(1-tail_fr,max_prob,n_point);

plot(icdf(tails_tot,plow_tail_tot),plow_tail_tot,'red','LineWidth',4);
plot(icdf(tails_tot,pint_tot),pint_tot,'black','LineWidth',2);
plot(icdf(tails_tot,pup_tail_tot),pup_tail_tot,'blue','LineWidth',2);

xlabel('Standardized Residuals')
ylabel('Probability')
title('Empirical Distribution Function TOT')
legend({'Pareto Lower Tail','Kernel Smoothed Interior','Pareto Upper Tail'},'location','northwest')
hold off

% ENI
subplot(2,2,4)
hold on 
grid on

min_prob = cdf(tails_eni,min(oil_std_returns.ENI_res));
max_prob = cdf(tails_eni,max(oil_std_returns.ENI_res));

plow_tail_eni = linspace(min_prob,tail_fr,n_point);
pint_eni = linspace(tail_fr,1-tail_fr,n_point);
pup_tail_eni = linspace(1-tail_fr,max_prob,n_point);

plot(icdf(tails_eni,plow_tail_eni),plow_tail_eni,'red','LineWidth',4);
plot(icdf(tails_eni,pint_eni),pint_eni,'black','LineWidth',2);
plot(icdf(tails_eni,pup_tail_eni),pup_tail_eni,'blue','LineWidth',2);

xlabel('Standardized Residuals')
ylabel('Probability')
title('Empirical Distribution Function ENI')
legend({'Pareto Lower Tail','Kernel Smoothed Interior','Pareto Upper Tail'},'location','northwest')
hold off

%% Asses the GPD Fit in the Lower Tails
figure('Name','GPD fit of the oil companies')
% CVX 
subplot(2,2,1)
[P,Q] = boundary(tails_cvx);
y = sortrows(oil_std_returns.CVX_res(oil_std_returns.CVX_res<Q(1),1))-Q(1);
plot(y,(cdf(tails_cvx,y+Q(1)))/P(1));
[F,x] = ecdf(y,'Bounds','on');
hold on
stairs(x,F,'r')
grid on

legend('Fitted GPD','Empirical CDF','Location','best')
xlabel('Exceedance')
ylabel('Probability')
title('Lower Tail of Standardized residuals of CVX')
ylim([0 1])

% COP 
subplot(2,2,2)
[P,Q] = boundary(tails_cop);
y = sortrows(oil_std_returns.COP_res(oil_std_returns.COP_res<Q(1),1))-Q(1);
plot(y,(cdf(tails_cop,y+Q(1)))/P(1));
[F,x] = ecdf(y);
hold on
stairs(x,F,'r')
grid on

legend('Fitted GPD','Empirical CDF','Location','best')
xlabel('Exceedance')
ylabel('Probability')
title('Lower Tail of Standardized residuals of COP')
ylim([0 1])

% TOT 
subplot(2,2,3)
[P,Q] = boundary(tails_tot);
y = sortrows(oil_std_returns.TOT_res(oil_std_returns.TOT_res<Q(1),1))-Q(1);
plot(y,(cdf(tails_tot,y+Q(1)))/P(1));
[F,x] = ecdf(y);
hold on
stairs(x,F,'r')
grid on

legend('Fitted GPD','Empirical CDF','Location','best')
xlabel('Exceedance')
ylabel('Probability')
title('Lower Tail of Standardized residuals of TOT')
ylim([0 1])

% ENI
subplot(2,2,4)
[P,Q] = boundary(tails_eni);
y = sortrows(oil_std_returns.ENI_res(oil_std_returns.ENI_res<Q(1),1))-Q(1);
plot(y,(cdf(tails_eni,y+Q(1)))/P(1));
[F,x] = ecdf(y);
hold on
stairs(x,F,'r')
grid on

legend('Fitted GPD','Empirical CDF','Location','best')
xlabel('Exceedance')
ylabel('Probability')
title('Lower Tail of Standardized residuals of ENI')
ylim([0 1])

%% Calculating the VAR with movable windows in order to compare and evaluate the method against the others.

TestWindowStart = find(year(oil_date.Date)==2020,1);
TestWindowEnd = find(year(oil_date.Date)==2022,1,'last');
TestWindow = TestWindowStart:TestWindowEnd;
WindowSize = 500;
index = oil_std_returns.Properties.VariableNames;
pVaR = [0.05,0.01];
VaREVT95 = zeros(length(TestWindow),width(index));
VaREVT99 = zeros(length(TestWindow),width(index));


%% Constructing the cycle for
for j = 1 : width(index)
    for t = TestWindow
        i = t-TestWindowStart+1;
        EstimationWindow = t-WindowSize:t-1;
        tails = paretotails(oil_std_returns{EstimationWindow,index(j)},tail_fr,1-tail_fr,'kernel');
        Q = boundary(tails);
        u = Q(1);
        beta = tails.LowerParameters(2);
        shape = tails.LowerParameters(1);
        N = WindowSize;
        Nu = sum(oil_std_returns{EstimationWindow,index(j)}<Q(1));
        VaREVT95(i,j) = EVT(u,beta,shape,N,Nu,pVaR(1));
        VaREVT99(i,j) = EVT(u,beta,shape,N,Nu,pVaR(2));
    end
    figure('Name','VaR EVT method');
    plot(oil_date.Date(TestWindow,1),[oil_std_returns{TestWindow,index(j)}]);
    hold on
    grid on
    plot(oil_date.Date(TestWindow,1),[-VaREVT95(:,j),-VaREVT99(:,j)]);
    xlabel('Date');
    ylabel('VaR EVT ');
    legend({'Returns','EVT 95','EVT 99'},'location','best');
    title(strcat(['VaR calculated with EVT for '],index(j)));
    % backtesting
    Returns_test = oil_std_returns{TestWindow,index(j)};
    Date_test = oil_date.Date(TestWindow,1);  
    vbt =varbacktest(Returns_test,[VaREVT95(:,j),VaREVT99(:,j)],'PortfolioID',index(j),'VaRID',{'EVT 95','EVT 99'},'VaRLevel',1-pVaR(1));
    backtest = runtests(vbt);
    disp(backtest)
end

%% Conclusion EVT
% PortfolioID     VaRID      VaRLevel      TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    ________    ________    ______    ______    ______    ______    ______    ______    ______    ______

    % "CVX_res"     "EVT 95"      0.95      yellow    reject    reject    accept    reject    accept    reject    reject
    % "CVX_res"     "EVT 99"      0.95      green     reject    reject    accept    reject    accept    reject    reject

% PortfolioID     VaRID      VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
 %   ___________    ________    ________    _____    ______    ______    ______    ______    ______    ______    ______

   %  "COP_res"     "EVT 95"      0.95      red      reject    reject    accept    reject    reject    reject    reject
    % "COP_res"     "EVT 99"      0.95      green    reject    reject    reject    reject    reject    reject    reject

% PortfolioID     VaRID      VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    ________    ________    _____    ______    ______    ______    ______    ______    ______    ______

    % "TOT_res"     "EVT 95"      0.95      red      reject    reject    accept    reject    reject    reject    reject
    % "TOT_res"     "EVT 99"      0.95      green    accept    accept    accept    accept    reject    reject    reject
        
% PortfolioID     VaRID      VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    ________    ________    _____    ______    ______    ______    ______    ______    ______    ______

     %"ENI_res"     "EVT 95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
     %"ENI_res"     "EVT 99"      0.95      green    reject    reject    accept    reject    accept    reject    reject

%% Our solution WS = 500
%PortfolioID     VaRID      VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
    %___________    ________    ________    _____    ______    ______    ______    ______    ______    ______    ______

    % "CVX_res"     "EVT 95"      0.95      green    accept    accept    accept    accept    accept    accept    accept
    % "CVX_res"     "EVT 99"      0.95      green    reject    reject    accept    reject    accept    reject    reject

    % PortfolioID     VaRID      VaRLevel      TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    ________    ________    ______    ______    ______    ______    ______    ______    ______    ______

     %"COP_res"     "EVT 95"      0.95      yellow    reject    reject    accept    reject    accept    reject    reject
     %"COP_res"     "EVT 99"      0.95      green     reject    reject    accept    reject    accept    reject    reject
    %  PortfolioID     VaRID      VaRLevel      TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
  %  ___________    ________    ________    ______    ______    ______    ______    ______    ______    ______    ______

   %  "TOT_res"     "EVT 95"      0.95      yellow    reject    reject    accept    accept    accept    reject    reject
    % "TOT_res"     "EVT 99"      0.95      green     reject    reject    accept    reject    accept    reject    reject
%  PortfolioID     VaRID      VaRLevel     TL       Bin       POF       TUFF       CC       CCI       TBF       TBFI 
   % ___________    ________    ________    _____    ______    ______    ______    ______    ______    ______    ______

    % "ENI_res"     "EVT 95"      0.95      green    accept    accept    accept    accept    reject    reject    reject
    % "ENI_res"     "EVT 99"      0.95      green    reject    reject    accept    reject    accept    reject    reject








