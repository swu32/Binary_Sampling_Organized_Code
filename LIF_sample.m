% input:
% D: conefficient for each neuron
% E: constant for each neuron, 1* n_neuron vector
% x1: lower bound for each neuron
% x2: upper bound for each neuron
%

% n_sample: numer of sampling steps
% bin: number of bins per sampling step
% Var: variance of gaussian noise
% V_b: baseline threshold
% samplelength: in ms, the number of sampling done
% dt: in unit of ms, timebin of each sampling unit

function [SpikeTrain,Ztrain,V_record,current_input,RX, P,samplelength, dt] = LIF_sample(D,E,x1,x2,R,A,C,n_sample,bin,Var,V_b,V_t,tau,method,P_pri,sigma,Y,G)
%% integrate and fire sampling process\
% refrac = 'on';% turning on refractory period

% in each sampling step it should be the same probability as the sampling
% process
time = n_sample; % ms
% bin = 100;
timebin = time*bin; % total number of bins in the loop
samplelength  = time;
dt = 1/timebin;
trial = 1;
% Var = 1;
Varbin = Var/bin;% variance per bin % for adding gaussian noise


% V_b = 0; % Baseline Voltage of -65mV
% V_t = 10; % Fireing threshold
tau = tau*bin; % membrane decay constant for 2 ms
V = V_b*ones(2,1); % give it a baseline voltage
% delta_t = 1/bin;
sdim = size(R,1);




SpikeTrain = zeros(sdim, trial,timebin);% matrix of spike train
Ztrain = zeros(sdim,trial,timebin); % EPSP current record
V_record = zeros(sdim,trial, timebin); % voltage record
u = ones(1,sdim);%?

% stepWidth = round(bin*3.06); % length of stepWidth
stepWidth = bin; % stepwidth should be the same length as sampling unit, same length as bin
refbin = 1;
% refbin = bin;
% beta = zeros(sdim,1);
% EPSP = [0,1;1,0];

EPSP = zeros(sdim,sdim);% - diag(ones(1,sdim));

% calculate b value when no current from other neurons
if method ==2
  beta = zeros(1,sdim);
  YTG = Y'*G;
  for k = 1:sdim % loop through neuron
    beta(k) = YTG(k)/(sigma^2)+R(k,k)/(2*sigma^2)+log(P_pri);
    %     EPSP(k,i)= get_EPSP(I, RF, k, i, sdim,S_pri,sigma,V_t, V_b, tau, bin,beta(k));
    EPSP(k,:)=R(k,:)./(sigma^2);
  end
  
  
elseif method == 1
  beta = C + A.*E;% beta = 1 x sdim
  % n^2 values to compute
  for k = 1:sdim % loop through neuron
    %     EPSP(k,i)= get_EPSP(I, RF, k, i, sdim,S_pri,sigma,V_t, V_b, tau, bin,beta(k));
    EPSP(k,:)=A.*D(k).*R(k,:);
  end
  
  
else % calculate b_k and EPSP ki using linear approximating current method
  YTG = Y'*G;
  beta = zeros(1,sdim);
  V_theta = V_t - V_b;

  for k = 1:sdim % loop through neuron
    r = P_pri/(1  - P_pri)*exp(YTG(k)/(sigma^2)+R(k,k)/(2*sigma^2));
    beta(k) = V_theta/(1 - exp((-1/r - 1)/tau));
    for q = 1:sdim
      EPSP(k,q) = V_theta*(exp((-1-r)/(r*tau)))/((1 - exp((-1-r)/(r*tau)))^2*r*tau)*R(k,q)/(sigma^2);      
    end
    
  end

    
  
end



% current_input = [n_neuron x total_bins]
% RX = R*states
% P: Probability at each bin
% V = [1 x n_neuron] temperorily stores voltage of previous values
% Ztrain = [n_neuron x total_trial x total_bins] stores the EPSP value of
% other neurons
%           temp1
YTG = Y'*G;
for j = 1:trial % sampling trials to run
  current_input = zeros(sdim,timebin);
  RX = zeros(sdim,timebin);
  P = zeros(sdim,timebin);
  
  
  for k = 1:sdim % initialization
    V(k) = V_b; % set to baseline threshold
  end
  
  for i = 1:timebin
    % input current for neuron1
    for k = 1:sdim
      
      Refractory = false;
      if i>=(refbin+1)
        if (sum(SpikeTrain(k,j,(i-refbin):i))>0)
          Refractory = true;
          V(k) = V_b; % set to baseline threshold
          V_record(k,j,i) = V(k);
          P(k,i) = 0;
          current_input(k,i) = 0;
        end
      else
        if (sum(SpikeTrain(k,j,1:i))>0)
          Refractory = true;
          V(k) = V_b; % set to baseline threshold
          V_record(k,j,i) = V(k);
          P(k,i) = 0;
          current_input(k,i) = 0;
        end
        
      end
      if Refractory == false % the neuron is not in refractory period  
         z = Ztrain(:,j,i);
          z = squeeze(z);

          temp1 = [z(1:k-1);1;z(k+1:end)];
          temp0 = [z(1:k-1);0;z(k+1:end)];
          value = R(k,:)*temp0;
          RX(k,i) = value;
        if method ==2
%         XTR_1 = temp1*R;
          P_unnorm_1 = exp(-1/(2*sigma^2)*( - 2*(YTG(k)) - 2.*R(k,:)*temp0 -R(k,k)))*(P_pri^1);
%           P_unnorm_1 = exp(EPSP(k,:)*temp1+beta(k) - log(A));
          P_unnorm_0 = (1 - P_pri);
          K = 1/(P_unnorm_1+P_unnorm_0);
%           P_t = K*P_unnorm_1;
          P_t = exp(beta(k)+EPSP(k,:)*temp0 + log(K));
          
          % calculate normalizing constant
          % (exp(E(k) - log(A)+D(k)*value) is ususally very small
          % verify current actually produce that probability
          u(k) = A*P_t+C;
          I = exp(bin./(P_t.*tau)).*(V_b - V_t)./(1-exp(bin./(P_t.*tau)));

          % need to check if log(K) is actually small
           % the probability that it is trying to achieve
          P(k,i) = P_t;
%           K
        elseif method == 1 
        % test whether in between x1(k) and x2(k)
        % obtain state of other neurons
        % \sum Rik zi
          
        if value <= x1(k)
          u(k) = C;
          P(k,i) = 0;
        elseif value >= x2(k)
          u(k) = A+C;
          P(k,i) = 1;
          
        else
          P_unnorm_1 = exp(-1/(2*sigma^2)*( - 2*(YTG(k)) - 2.*R(k,:)*temp0 -R(k,k)))*(P_pri^1);
          P_unnorm_0 = (1 - P_pri);
          K = 1/(P_unnorm_1+P_unnorm_0);
          P_t = K*P_unnorm_1;
          u(k) = beta(k) + EPSP(k,:)*temp0; % add current from other neurons
          P_linear = D(k)*value+E(k);
          
          
          P(k,i) = D(k)*value+E(k);
        end
        
        else
          u(k) = beta(k) + EPSP(k,:)*temp0; % add current from other neurons
          P_unnorm_1 = exp(-1/(2*sigma^2)*( - 2*(YTG(k)) - 2.*R(k,:)*temp0 -R(k,k)))*(P_pri^1);
          P_unnorm_0 = (1 - P_pri);
          K = 1/(P_unnorm_1+P_unnorm_0);    
          P_t = K*P_unnorm_1;
%           display(P(k,i))
          % compare for a sanity check
          P(k,i) = (tau*log(u(k)/(u(k)+(V_b - V_t))))^(-1);
        end
        
        
        
        current_input(k,i) = u(k);
        
        
        V(k) = V_b + u(k) + exp(-1/tau)*(V(k)- V_b - u(k));
        %       V(k) = V(k)+ normrnd(0,sqrt(Varbin));
        
        if V(k)<V_b
          V(k) = V_b;
        end
        
        V_record(k,j,i) = V(k);
        
        
        
        if V(k)>V_t
          SpikeTrain(k,j,i) = 1;
          Ztrain(k,j,i) = 1;
          for l = 1:stepWidth
            if (i+l)>timebin
              Ztrain(k,j,i:end) = 1;
            else
              Ztrain(k,j,i+l) =1;
            end
          end
          V(k) = V_b;% Set Membrane Potential Back to baseline after firing
          V_record(k,j,i) = 20;
          
        end
        
      end
      
      
    end
  end
end
