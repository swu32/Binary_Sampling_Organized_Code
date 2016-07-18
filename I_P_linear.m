% SHUCHEN WU
% June 2016
% linear approximation of I as a function of p
% I = a*P_k+c;


% linear fit doesnt work very well when P < 0.1
% should be an optional statement when p<0.1
% obtain I's coefficient

% input: 

% I_P_linear(tau, V_r, V_t,plot)


function [A,C] = I_P_linear(varargin)
bin = 1000;
if nargin > 1
  tau = varargin{1}*bin;
  V_r = varargin{2};
  V_theta = varargin{3};
  makeplot = varargin{4}; 
  intersect_zero = varargin{5};

else
makeplot = false;
V_r = 0;
V_theta = 15;
tau = 10*bin;
end
P = (0.01:0.01:1);
I = exp(bin./(P.*tau)).*(V_r - V_theta)./(1-exp(bin./(P.*tau)));
P_fit = polyfit(P,I,1);% linear fit I as a linear function of P

A = P_fit(1);
C = P_fit(2);

if intersect_zero
  A = I(:)\P(:);
  C = 0;
  
end



if makeplot
  subplot(2,1,1)
  yfit = P_fit(1)*P+P_fit(2);
  plot(P,I,'o')
  hold on
  plot(P,yfit,'r-');
  xlabel('P')
  ylabel('I')
  legend('theortical','linear fit')
  title('fit current as a linear function of probability')
  
  subplot(2,1,2)
  residual = I - yfit;
  plot(P,residual,'o')
  xlabel('probability')
  ylabel('residual, I - fitline')
  title('Linear Fit Residual')

  
  figure
  yfit = P_fit(1)*P+P_fit(2);
  plot(I,P,'o')
  hold on
  plot(yfit,P,'r-');
  xlabel('I')
  ylabel('P')
  legend('theortical','linear fit')
  title('Inverse Fit')

  
end



