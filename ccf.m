function ta = ccf(x,y,p)
% CCF - Compute Correlogram Through p Lags

% --------------------------
% USER INPUT CHECKS
% --------------------------

[n1, n2] = size(y) ;
if n2 ~=1
    error('Input series y must be an nx1 column vector')
end

[a1, a2] = size(p) ;
if ~((a1==1 & a2==1) & (p<n1))
    error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
end



% -------------
% BEGIN CODE
% -------------

ta = zeros(p,1) ;
global N 
N = max(size(y)) ;
global ybar 
global xbar
ybar = mean(y); 
xbar = mean(x);


% Collect ACFs at each lag i
for i = 1:p
   ta(i) = acf_k(x,y,i) ; 
end

% Plot ACF
% Plot rejection region lines for test of individual autocorrelations
% H_0: rho(tau) = 0 at alpha=.05
bar(ta)
line([0 p+.5], (1.96)*(1/sqrt(N))*ones(1,2))
line([0 p+.5], (-1.96)*(1/sqrt(N))*ones(1,2))

% Some figure properties
line_hi = (1.96)*(1/sqrt(N))+.05;
line_lo = -(1.96)*(1/sqrt(N))-.05;
bar_hi = max(ta)+.05 ;
bar_lo = -max(ta)-.05 ;

if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
    axis([0 p+.60 line_lo line_hi])
else
    axis([0 p+.60 bar_lo bar_hi])
end
title({' ','Correlogram',' '})
xlabel('Lag Length')
set(gca,'YTick',[-1:.20:1])
% set number of lag labels shown
if (p<28 & p>4)
    set(gca,'XTick',floor(linspace(1,p,4)))
elseif (p>=28)
    set(gca,'XTick',floor(linspace(1,p,8)))
end
set(gca,'TickLength',[0 0])




% ---------------
% SUB FUNCTION
% ---------------
function ta2 = acf_k(x,y,k)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% 
global ybar
global xbar
global N
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(x(i-k)-xbar) ;
end

% Denominator, unscaled variance
var =sqrt(sum((y-ybar).^2)*sum(((x-xbar).^2))) ;

ta2 = sum(cross_sum) / var ;

