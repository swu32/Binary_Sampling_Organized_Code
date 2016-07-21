% gives two bounds for curve fitting in section
% Fit the curve using piecewise function, below x1 and above x2 are 0 or 1
% P = D*(sum( RFki*Zi )) + E;% linear as a funtion of other neurons
% y = P
% x = (sum( RFki*Zi ))

% input:   array to fit
% output:  D   linear coefficient
% E:  coefficient constant
% bd: a bound to sort out the rest, if too tight, no points will get into the middle

% 
function [D,E,x1,x2,x_fit,y_fit] = cf_fit(x,y,bd,PLOT)
% there is a another possbility of being Nan
Na = isnan(y);
x = x(~Na);
y = y(~Na);


x = sort(x);
y = sort(y);
A0 = y<=bd;% close to 0
A1 = y >= (1.00-bd);% close to 1
% A = (y>=(0+bd)) & (y <= (1-bd));
A = ~(A0+A1);
SUM = sum(A) + sum(A0) +sum(A1)
sum(A)
sum(A0)
sum(A1)


% find D,E, best fits 
if sum(A)>0
  % fit those parts and find x1,x2,D,E
% obtain linear fitting coefficients
  [fitvars, resNorm]= polyfit(x(A), y(A), 1);

  D = fitvars(1);
  E = fitvars(2);


  % range of linear fitting
  x1 = min(x(A));
  x2 = max(x(A));

  x_f = x(A);
  y_f = polyval(fitvars,x_f);
  % for those smaller than 0 or larger than 1, fix them to 0 or 1
  B_0 = y_f<0;
  y_f(B_0) = 0;
  B_1 = y_f>1;
  y_f(B_1) = 1;
  
  x_f0 = x(A0);
  y_f0 = 0.*x_f0;
  
  x_f1 = x(A1);
  y_f1 = 1.*ones(1,length(x_f1));
  
  
  y_fit = zeros(1,length(y));
  y_fit(A) = y_f;
  y_fit(A0) = y_f0;
  y_fit(A1) = y_f1;
  
  x_fit = zeros(1,length(x));
  x_fit(A) = x_f;
  x_fit(A0) = x_f0;
  x_fit(A1) = x_f1;
  
  
  
elseif sum(A)==0 % if all zeros or all ones or either zero or 1, nothing in between, set D and E to 0.   
  sum(A)
      D = 0;
      E = 0;
      
  if sum(A0)==length(x)% all zero
    x_fit = x;
    y_fit = zeros(1,length(y));

    x1 = 1000;
    x2 = 2000;
    
  elseif sum(A1)==length(x)
    x2 = -1000;
    x1 = -2000;
    x_fit = x;
    y_fit = ones(1,length(y));
    
  else% either 0 or 1. 
    x_fit = zeros(1,length(x));
    x1 = x((A0));
    x1 = x1(end);
    x2 = x((A1)); 
    x2 = x2(1);
    x_f0 = x(A0);
    y_f0 = 0.*x_f0;
    
    
    x_f1 = x(A1);
    y_f1 = 1.*ones(1,length(x_f1));

    x_fit(A0) = x_f0;
    x_fit(A1) = x_f1;
    y_fit(A0) = y_f0;
    y_fit(A1) = y_f1;
    
    
  end
  
  
end



% TF = isempty(x1);
% if TF % check whether it is all one, all zero, or either one or zero
% allzero = true;
% allone = true;
% K = y;
% 
% % maybe filter out negative slope 
% for i = 1:length(K)
%   if K(i)>bd & K(i)<(1-bd)
%     allzero = false;
%     allone = false;
%   elseif K(i)<BD
%     allone = false;
%   elseif K(i)>(1-BD)
%     allzero = false; 
%   end
% end
% 
% if allzero % set x1 to infinity
%   x1 = 1000;
%   x2 = 2000;
% end
% 
% if allone
%   x2 = -1000;
%   x1 = -2000;
% end  
%   
%   if (~allzero) & (~allone) % both false, but nothing falls in inteval
%     % take x1 and x2 as in the intersection between 1 and 0
%     x1 = X(min(find(Y>0.5)));
%     x2 = X(max(find(Y<0.5)));
%   end
%   
%   
%   
% else
% % fit the interval between x1 and x2
% %  x_f = x1:0.01:x2;
% x_f = x(A);
% y_f = polyval(fitvars,x_f);
% % for those smaller than 0 or larger than 1, fix them to 0 or 1
% B_0 = y_f<0;
% y_f(B_0) = 0;
% B_1 = y_f>1;
% y_f(B_1) = 1;
% 
% x_f0 = x(A0);
% y_f0 = 0.*x_f0;
% 
% x_f1 = x(A1);
% y_f1 = 1.*ones(1,length(x_f1));
% 
% 
% y_fit = zeros(1,length(y));
% y_fit(A) = y_f;
% y_fit(A0) = y_f0;
% y_fit(A1) = y_f1;
% 
% x_fit = zeros(1,length(x));
% x_fit(A) = x_f;
% x_fit(A0) = x_f0;
% x_fit(A1) = x_f1; 
%   
%   
% end





x1
x2
D
E
if strcmp(PLOT,'on')
figure


%
subplot(211)
% 
plot(x,y,'o',x_fit,y_fit,'+')
title('Plot of Data (Points) and Model (Line)')

% plot residual
res = y - y_fit;

subplot(212)
plot(x,res,'+')
title('Plot of the Residuals')
end

  