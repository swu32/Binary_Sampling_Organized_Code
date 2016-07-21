% SHUCHEN WU
% June 2016
% main program to run Sampling in comparison with LIP neurons

% method 1: piecewise linear fit, linear readout
% method 2: exponential fit
% method 3: direct taylor expansion of current

figure
close all;
fclose('all');
% exp_I = 'off';
method = 1;
% obtain a linear approximation of each neuron to its response

% it is important to obtain natural image and receptive fields, use the
% latest learned 128 receptive fields


% do one image as first step

% obtain linear approximation of function close to mean

% record for every single neuron and run for LIF, compare result with
% samplng


%% load a set of basis functions
[G, P_pri,sigma] = ReadDataFile; 
pixel = size(G,1);
xdim = pixel;

% extract an image from von hatern dataset
Y = extractImagePatches( 1, round(sqrt(pixel)) );

% I = rand(8,8);



%% fix coefficient and obtain coefficient
tau = 1;% ms
V_r = 0;
V_theta = 15;% threshold voltage
bin = 100;
Var = 0;% variance of gaussian noise if it is being added. 

% I_P_linear(tau, V_r, V_t,plot)


%TO DO: Evaluate definite integral solution and compare with the matlab fit
[A,C] = I_P_linear(tau, V_r, V_theta,false,false);  % obtain coefficient of I as a function of P, only needs to be computed once

% set up the number of neurons
sdim = 10;
GG = G(:,1:sdim);
% R = -GG'*GG;
R = getR(GG);
YTG = (Y'*GG);



%% enumerate all possible combinatorics of a Binary Block code with length
% equal to the number of neurons sdim
% instead of enumberating all possible combinatorics, generate random
% numbers instead. Will approixmate given law of large numbers

% number of random numbers generated
n_rand = 100000;


x = zeros(sdim,n_rand); % record every x value, conditional probability when x is on
y = zeros(sdim,n_rand); % record exp(x something)
D = zeros(1,sdim);% first degree coeffienct
E = zeros(1,sdim);% constant

if method == 2


  % use exponential function to express current
  % obtain theortical calculation of D and E.
  %   D = (1/(2*sigma^2)).*ones(1,sdim);
for k  = 1:sdim
  E(k) =  YTG(k)/(sigma^2) +R(k,k)/(2*sigma^2)+ log(P_pri);
end
%   E = log(A) + (1/sigma^2).*(Y'*GG) + log(P_pri);
  D = (1/(sigma^2)).*ones(1,sdim);
  x1 = [];
  x2 = [];
  
elseif method ==1
  
x1 = zeros(1,sdim);% lower bound for linear fit
x2 = zeros(1,sdim);% upper bound 
resNorm = zeros(1,sdim);

% loop through every single neuron, get every single conditional
% probability % better not to do that

% TODO: check why it is different from sampling result!!!
% for each neuron obtain a linear approximation


for i = 1:sdim 
  for j = 1:n_rand % loop through all states of other neurons
%     B = rand(1,sdim)>(1-1/sdim); % with a small probability of firing
    B = round(rand(1,sdim));
    prevs = double((B));
    
    % check probability calculated from derived equation
    temp1 = [prevs(1:i-1),1,prevs(i+1:end)];
    temp0 = [prevs(1:i-1),0,prevs(i+1:end)];
    % i is the dimension to be 
    L2_s_1 = sum((Y(:)-GG*temp1').^2);% 64x1 (G*temp1')^2 = (G*temp1')'*(G*temp1') = temp1*G'*G*temp1'
    L2_s_0 = sum((Y(:)-GG*temp0').^2);
    
    post_s_1 = exp(-1/(2*sigma^2)*(L2_s_1))*(P_pri^1);
    post_s_0 = exp(-1/(2*sigma^2)*(L2_s_0))*(1-P_pri)^1;
    C = 1/(post_s_1+post_s_0);
    
    P_post_s1 = C*post_s_1;
    P_post_s0 = C*post_s_0;
    % case both values are close to 0, 
    if isinf(C)
      if post_s_1 > post_s_0
        P_post_s1 = 1;
        P_post_s0 = 0;
      else
        P_post_s0 = 1;
        P_post_s1 = 0;
      end      
    end
        
    % check probability calculated from derived equation
    
    post_unnorm_1 = exp(-1/(2*sigma^2)*( - 2*(YTG(i)) - 2.*R(i,:)*temp0' - R(i,i)))*(P_pri^1);
    post_unnorm_0 = ((1 - P_pri)^1);        

    K = 1/(post_unnorm_1+post_unnorm_0);
    post_1 = K*post_unnorm_1;
    
    
    % check probability calculated from ratio
    r = P_pri/(1  - P_pri)*exp(YTG(i)/(sigma^2)+R(i,i)/(2*sigma^2) + R(i,:)*temp0'/(sigma^2));
    P = r/(r+1);
    
    if abs(P_post_s1-post_1)>10e-04
           K

    end
        
    
    y(i,j) = P_post_s1;    
    x(i,j) = R(i,:)*temp0';% sum(R_jk*x_j)
    % calculate variable as a function of feedforward input and
    % recurrent connection
            
    
  end
%   Na = isnan(y(i,:));
 eps = 10e-3;
 [D(i),E(i),x1(i),x2(i),x_fit,y_fit] = cf_fit( x(i,:),y(i,:),eps,'off');
 % sanity check comparison:
 Y1 = sigma^2*log(eps/(1-eps)*((1-P_pri)/P_pri)) - R(i,i)/2-YTG(i);
 Y2 = sigma^2*log((1-eps)/eps*((1-P_pri)/P_pri)) - R(i,i)/2-YTG(i);
 YY = linspace(Y1,Y2,1000);
 % obtain approximation of int_I and int_pI
 r = P_pri/(1-P_pri)*exp(YTG(i)/(sigma^2) + YY./(sigma^2) + R(i,i)/(2*sigma^2));
 P = 1./(1+1./r);
 value_int_P = P;
% value_int_I = 3.*p+2;

value_int_P = value_int_P(2:end);
value_int_P = value_int_P.*((Y2 - Y1)/(length(YY)));

int_P = sum(value_int_P(~isnan(value_int_P)));


% value_int_pI = (3.*p.^2+2.*p);
value_int_YP = YY.*P;
value_int_YP = value_int_YP(2:end);
value_int_YP = value_int_YP.*((Y2 - Y1)/(length(YY)));
int_YP = sum(value_int_YP(~isnan(value_int_YP)));

% left * M = right
% right = [int_I  ]
%         [int_pI ];
% M = [ a  b ]
%     [ c  d ];
% left = [ C ]
%        [ A ];

right = [int_P;int_YP];
a = Y2 - Y1;
b = 0.5*(Y2^2 - Y1^2);
c = 0.5*(Y2^2 - Y1^2);
d = 1/3*(Y2^3 - Y1^3);
M = [a,b;c,d];
left = (M)\right;

EE = left(1);
DD = left(2);

% caution the possibility of resulting in a negative probability
% sanity check
figure

DEfit  = zeros(1,n_rand);

for k = 1:n_rand
  if (x(i,k)<Y1)
    DEfit(k) = 0;
  elseif(x(i,k)>Y2)
    DEfit(k) = 1;
  else
    DEfit(k) = DD*x(i,k)+EE;
    if DEfit(k)<0
      DEfit(k) = 0;
    elseif DEfit(k)>1
      DEfit(k) = 1;        
    end
    
  end
end



%
subplot(211)
% 
plot(x(i,:),y(i,:),'o',x_fit,y_fit,'+')
hold on
plot(x(i,:), DEfit,'+')
axis([Y1 Y2 -0.2 1.2])
title('Plot of Data (Points) and Model (Line)')
% res2 = [zeros(): :ones()];
% plot residual
res1 = y(i,:) - y_fit;
res2 = y(i,:) - DEfit;
subplot(212)
plot(x(i,:),res1,'+')
hold on
plot(x(i,:),res2,'+')

title('Plot of the Residuals')
 

end


  % obtain curive fitting range
  else method == 3

    X = rand(n_rand,sdim)>(1-1/sdim); % with a small probability of firing
    X = double((X));
    
    
    
for i  = 1:sdim
  for j = 1:size(X,1) % loop through all states of other neurons
    x_on = [X(j,1:i-1),1,X(j,i+1:end)];
    x_off = [X(j,1:i-1),0,X(j,i+1:end)];
    % For the sake of comparison and sanity check
    L2_s_1 = sum((Y-GG*x_on').^2);% For a not so good image, L2 can be huge
    L2_s_0 = sum((Y-GG*x_off').^2);
    post_s_1 = exp(-1/(2*sigma^2)*(L2_s_1))*(P_pri^1);
    post_s_0 = exp(-1/(2*sigma^2)*(L2_s_0))*(1-P_pri)^1;
    P_post_s1 = 1/(post_s_1+post_s_0)*post_s_1;
    
    
    % 
    P_2_1 = P_pri*exp(YTG(k)/(sigma^2)+R(k,k)/(2*sigma^2) + R(k,:)*x_off'/(sigma^2));
    P_2_0 = 1-P_pri;
    K =1/(P_2_1+P_2_0);
    P_2 = K*P_2_1;
    
    %          
    post_unnorm_1 = exp(-1/(2*sigma^2)*( - 2*(YTG(k)) - 2.*R(k,:)*x_off' -R(k,k)))*(P_pri^1);
    post_unnorm_0 = exp(-1/(2*sigma^2)*(0))*((1 - P_pri)^1);   
    KK = 1/(P_unnorm_1+P_unnorm_0);
    P_t = KK*P_unnorm_1;
    
    % real shit
    r = P_pri/(1  - P_pri)*exp(YTG(k)/(sigma^2)+R(k,k)/(2*sigma^2) + R(k,:)*x_off'/(sigma^2));
    P = r/(r+1);
    I_k = (V_theta)/(1-exp((-1-r)/(tau*r)));
    
    y(i,j) = I_k;
    x(i,j) = R(i,:)*x_off';% sum(R_jk*x_j)    


  end
    
end    
    
end


% do many times of this process to obtain general statistics
n_run = 1;
P_LIP = zeros(n_run,sdim);
P_Sample = zeros(n_run,sdim);
P_Sample_2 = zeros(n_run,sdim);

for i = 1:n_run
  
%% Do Sampling at this step
n_sample = 100;
% DO two times of sampling step for compariosn of its own statistics

[S_1, P_e_s_1,P_post_s1_1] = sampling(R,GG, P_pri,sigma,Y,n_sample,1,1,'off');
[S_2, P_e_s_2,P_post_s1_2] = sampling(R,GG, P_pri,sigma,Y,n_sample,1,1,'off');
S = S_1;
P_e_s = P_e_s_1;
P_post_s1 = P_post_s1_1;

disp('Hello World')


%% Do integrate and fire Process

[SpikeTrain,Ztrain,V_record,current_input,RX,P,samplelength, dt] = LIF_sample(D,E,x1,x2,R,A,C,n_sample,bin,Var,V_r,V_theta,tau,method,P_pri,sigma,Y,GG);



P_LIP(i,:) = squeeze(sum(SpikeTrain,3))./size(S,2);

P_Sample(i,:) = squeeze(sum(S,2))./size(S,2);
P_Sample_2(i,:) = squeeze(sum(S_2,2))./size(S_2,2);



end


%% Statistics
% check sampling probability
close all


% raster spike plot
figure
tVec = linspace(1,n_sample,n_sample);
plotRaster(S,tVec,'b')
tVec = linspace(1,n_sample,n_sample*bin);
plotRaster(SpikeTrain,tVec,'r')

xlabel('time')
ylabel('neuron')
title('100 neuron, tau = 1')
% sampling blue, LIP red

% raster spike plot next to each other
figure
subplot(1,2,1)
tVec = linspace(1,n_sample,n_sample);
plotRaster(S,tVec,'b')
title(strcat('sampling, 10 neurons, tau',num2str(tau)))
xlabel('time')
ylabel('neuron')

subplot(1,2,2)
tVec = linspace(1,n_sample,n_sample*bin);
plotRaster(SpikeTrain,tVec,'r')
title(strcat('LIF, 10 neurons, tau',num2str(tau)))
xlabel('time')
ylabel('neuron')
% sampling blue, LIP red


% compare posterior probability: 

%scatter plot
figure
% P_LIP = squeeze(sum(SpikeTrain,3))./size(SpikeTrain,3);
% P_LIP = squeeze(sum(SpikeTrain,3))./size(S,2);
% 
% P_Sample = squeeze(sum(S,2))./size(S,2);
% P_Sample_2 = squeeze(sum(S_2,2))./size(S_2,2);

x =  (0:0.001:1);
plot(x,x);
hold on
% plot(P_Sample, P_Sample_2,'o')
% hold on
% plot(P_Sample, P_LIP,'o')
plot(P_Sample(:), P_Sample_2(:),'o')
hold on
plot(P_Sample(:), P_LIP(:),'o')
xlabel('Sampling Probability')
ylabel('LIP Probability')
legend('ideal','sample1-sample2','sample-LIP')
title(strcat('Scatter Plot of Marginal Probability with tau = ',num2str(tau)))



% compare distribution from the y = x line eith sampling
% dist_sample = abs(P_Sample - P_Sample_2)./sqrt(2);
% dist_LIP = abs(P_Sample - P_LIP)./sqrt(2);
dist_sample = (P_Sample(:) - P_Sample_2(:))./sqrt(2);
dist_LIP = (P_Sample(:) - P_LIP(:))./sqrt(2);
figure
[x1,y1] = hist(dist_sample,10);
% bar(y1,x1,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
hold on
[x2,y2] = hist(dist_LIP,10);
% bar(y2,x2,'FaceColor',[.5 0 .5],'EdgeColor',[0 1 1],'LineWidth',1.5)
bpcombined = [x1(:),x2(:)];
hb = bar(y2, bpcombined);
legend('sampling-sampling distribution','LIP-sampling distribution')
xlabel('distance from ideal y = x line')
ylabel('count')

% Debugging statistics
figure
% S: Sampling Spike
% P_Post_s1: Posterior probability of Spiking
for i = 1:sdim
  subplot(sdim,1,i);
  plot(S(i,:),'o');
  hold on;
  plot(P_post_s1(i,:),'-o');
end

% compare sampling and LIP spike probability
figure
for i = 1:sdim
subplot(sdim,1,i);
t = linspace(0,n_sample,length(P(i,:)));
plot(t,P(i,:),'o');
hold on
t = linspace(0,n_sample,length(P_post_s1(i,:)));
plot(t,P_post_s1(i,:),'-o');
end



% plot input current
% probability of spiking. 
% x1, x2
% R_ik Xi
figure
for i = 1:sdim
  subplot(sdim,1,i);
  plot(x1(i).*ones(1,size(RX,2)),'b')
  hold on
  plot(x2(i).*ones(1,size(RX,2)),'b')
  hold on
  plot(RX(i,:),'r')
%   hold on
%   plot(squeeze(P(i,:)),'k')
  legend('x1','x2','Rki_zi')
end



figure
for i = 1:sdim
  subplot(sdim,1,i);
  plot((A+C)*ones(1,size(RX,2)),'b')
  hold on
  plot(C*ones(1,size(RX,2)),'b')
  hold on
  plot(current_input(i,:),'r')
  plot(squeeze(P(i,:)),'k')

  legend('A+C','A','current input','probability')
end

% plot SpikeTrain: LIF Spike
% V_record: LIF voltage record
% P: Probability of Spike in each time frame, every time neurons state
% change, P is going to change
figure
title('SpikeTrain and voltage record for LIP neuron sampling process')

for i = 1:sdim
  subplot(sdim,1,i);
  semilogy(squeeze(SpikeTrain(i,:)),'o');
  hold on
  semilogy(squeeze(V_record(i,1,:)))
%   hold on
%   plot(squeeze(P(i,:)),'k')
  legend('SpikeTrain','voltage record')
end


% check if reconstruced RF are close to original patterns
patch = 1;
figure
subplot(131)
imagesc(reshape(Y(:),sqrt(xdim),sqrt(xdim)))
title('original patch')
subplot(132)
imagesc(reshape(GG*SpikeTrain(:,patch,end),sqrt(xdim),sqrt(xdim)))
title('end of trial reconstruction')
subplot(133)
E = mean(GG*squeeze(SpikeTrain(:,patch,100:end)),2);
imagesc(reshape(E,sqrt(xdim),sqrt(xdim)))
title('averaged reconstruction after 1ms')
colormap('gray')


% reconstruction error comparing LIF and Sampling
subplot(1,2,1)
grad = 100;
Mean_err = zeros(1,(size(SpikeTrain,3)/grad));
Median_err = Mean_err;
sample = Median_err;
for i = 1:(size(SpikeTrain,3)/grad)
  % square reconstruction error
A = (abs(Y - GG*squeeze(SpikeTrain(:,1,i*grad))));
% A = (sum(A,1))./(sdim);
Mean_err(i) = mean(A);
Median_err(i)=median(A);
sample(i) = i;
end
plot(sample,Mean_err,'-o')
hold on
plot(sample,Median_err,'-o')
legend('Mean LIF Error','Median LIF Error')

title('Reconstruction Error across LIF Stages')
ylabel('reconstruction error per pixel')

xlabel('LIF step x100')

grad = 1;
subplot(1,2,2)
Mean_err_s = zeros(1,(size(S,2)/grad));
Median_err_s = zeros(1,(size(S,2)/grad));
sample_s = zeros(1,(size(S,2)/grad));
for i = 1:(size(S,2)/grad)
  % square reconstruction error
A = (abs(Y - GG*squeeze(S(:,i*grad))));
% A = (sum(A,1))./(sdim);
Mean_err_s(i) = mean(A);
Median_err_s(i)=median(A);
sample_s(i) = i;
end
plot(sample_s,Mean_err_s,'-o')
hold on
plot(sample_s,Median_err_s,'-o')
legend('Mean Sampling Error','Median Sampling Error')

title('Reconstruction Error across Sampling Stages')
ylabel('reconstruction error per pixel')

xlabel('sampling step x1')






figure
subplot(131)
imagesc(reshape(Y(:),sqrt(xdim),sqrt(xdim)))
title('original patch')
subplot(132)
imagesc(reshape(GG*S(:,end),sqrt(xdim),sqrt(xdim)))
title('end of trial reconstruction')
subplot(133)
E = mean(GG*squeeze(S(:,1:end)),2);
imagesc(reshape(E,sqrt(xdim),sqrt(xdim)))
title('averaged reconstruction after 1ms')
colormap('gray')

% S: Sampling Spike
% SpikeTrain: LIF Spike

isi_vect = [];
for i = 1:sdim
  Times = dt*find(abs(squeeze(SpikeTrain(i,1,:))'-1) < 0.00000001);
  isi_vect = [diff(Times),isi_vect];
end
% dt in unit of ms

% Distribution of ISI
figure
hist(isi_vect,200);
xlabel('time(ms)')
ylabel('count')
title(strcat('histogram of interspike interval distribution with tau = ',num2str(tau)))




figure
Mean = mean(SpikeTrain(:,1,:),3);
Var = var(SpikeTrain(:,1,:),0,3);
plot(Mean(:),Var(:),'+')
xlabel('Mean spikes per 1/100 ms')
ylabel('Spike Variance per 1/100 ms')
title(strcat('Mean and Variance of Spike Train for LIF neuron with tau = ',num2str(tau)))




figure
Fano = Var(:)./Mean(:);
hist(Fano(:))
title('histogram of Fano factor distribution per 1/100 ms')



disp('Hello World')


%plot conditional probability as a function of first two terms