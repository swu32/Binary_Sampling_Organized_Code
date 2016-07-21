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
% Sanity checked and it worked!

% P = (0.01:0.01:1);
% I = exp(bin./(P.*tau)).*(V_r - V_theta)./(1-exp(bin./(P.*tau)));
% P_fit = polyfit(P,I,1);% linear fit I as a linear function of P

% A = P_fit(1);
% C = P_fit(2);


% obtain approximation of int_I and int_pI

p = linspace(0,1,10000);
value_int_I = exp(bin./(p.*tau)).*(V_r - V_theta)./(1-exp(bin./(p.*tau)));
% value_int_I = 3.*p+2;

value_int_I = value_int_I(2:end);
value_int_I = value_int_I.*(1/(length(p)));

int_I = sum(value_int_I(~isnan(value_int_I)));


% value_int_pI = (3.*p.^2+2.*p);
value_int_pI = p.*exp(bin./(p.*tau)).*(V_r - V_theta)./(1-exp(bin./(p.*tau)));
value_int_pI = value_int_pI(2:end);
value_int_pI = value_int_pI.*(1/(length(p)));
int_pI = sum(value_int_pI(~isnan(value_int_pI)));

% left * M = right
% right = [int_I  ]
%         [int_pI ];
% M = [ a  b ]
%     [ c  d ];
% left = [ C ]
%        [ A ];

right = [int_I;int_pI];
a = 1;
b = 1/2;
c = 1/2;
d = 1/3;
M = [a,b;c,d];
left = (M)\right;

C = left(1);
A = left(2);




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



