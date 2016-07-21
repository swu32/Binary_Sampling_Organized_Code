close all;
%% Sampling
% Inference on state of S using Gibbs Sampling
% knowing S in a certain ground truth, the probability of S being in
% certain state

[W,P_pri,sigma] = ReadDataFile;
% sdim Dimension of coefficients
% xdim dimension for image(number of pixels)
[xdim,sdim] = size(W);

% W = randn(xdim,sdim);% Learned Weights


n_patches = 50;%
n_sample = 1000; % number of samples taken by s
% sigma = 0.1;

% S has dimension: n_coefficient*n_samples*n_patches
S = zeros(sdim,n_sample,n_patches);% one or 0
L2_1 = zeros(sdim, n_sample, n_patches);
L2_0 = zeros(sdim, n_sample, n_patches);

% X = W*S + sigma*randn(xdim,samples);

% TODO;
% loadsample turn out to be sqrt(xdim)*sqrt(xdim)*n_samples
% load learned weights

X = extractImagePatches( n_patches, sqrt(xdim) );


% P_pri = 0.3;% probability of prior to be 1
S_pri = ones(sdim,1)*P_pri;% probability for each neuron to be active

%posterior 
P_post_s_1 = ones(sdim,n_sample,n_patches);
P_post_s_0 = ones(sdim,n_sample,n_patches);
P_post_s_1(:,1,1) = S_pri;
P_post_s_0(:,1,1) = 1-S_pri;

% initialize S
Initialprevs = zeros(1,sdim);
for i = 1:sdim
  if rand(1)>P_pri
    Initialprevs(i) = 0;
    S(i,1,:) = 0;% update the first set of S for all patches  
  else
    Initialprevs(i) = 1;
    S(i,1,:) = 1;
  end
end


% Sample from image patches
for j = 1:n_patches% sees new image
  
  
  % Everytime seeing a new image, update the initialization of prevs
  prevs = Initialprevs;% update initialization of S when images is renewed
  
  for i = 2:n_sample   
    for k = 1:sdim
      % L2 norm
      % different values of L2 assuming Sk is 0 and 1
      L2_s_1 = sum((X(:,j)-W*[prevs(1:k-1),1,prevs(k+1:end)]').^2);% For a not so good image, L2 can be huge
      L2_s_0 = sum((X(:,j)-W*[prevs(1:k-1),0,prevs(k+1:end)]').^2);
      L2_1(k,i,j) = L2_s_1;
      L2_0(k,i,j) = L2_s_0;
      % for comparison:
      
      
      

      % Calculate Posterior
      % in Big L2, post_s_1 is almost 0;
      post_s_1 = exp(-1/(2*sigma^2)*(L2_s_1))*(S_pri(k)^1);
      post_s_0 = exp(-1/(2*sigma^2)*(L2_s_0))*(1-S_pri(k))^1;
      
      P_post_s_1(k,i,j) = 1/(post_s_1+post_s_0)*post_s_1;
      P_post_s_0(k,i,j) = 1/(post_s_1+post_s_0)*post_s_0;
      
      if(rand(1)> P_post_s_1(k,i,j))
        S(k,i,j) = 0;
        prevs(k) = 0;
      else
        S(k,i,j) = 1;
        prevs(k) = 1;
      end      
      % Normalize Posterior      
      % update the value for sprev      
    end
  end
end



% normalize by patch variance and add preprocessing step for inference and learning
% run a simulation in which you do all the plots for an undercomplete case  
% add cross-correlogram plots between the neurons who?s PFs agree the most with each other
% restrict auto-corr plot to relevant range.
% compute FF using 1-p, not var/mean, to account for p?s so small no 1?s are sampled


%% Analysis



% check if reconstruced RF are close to original patterns
figure(1)
n_patches = 1;
patch = 6;
% calculate expectation value
E = mean(W*S(:,100:end,patch),2);
subplot(n_patches,3,1)
imagesc(reshape(X(:,patch),sqrt(xdim),sqrt(xdim)))
% title('sample original patch')
title('a')
subplot(n_patches,3,2)
imagesc(reshape(W*S(:,999,patch),sqrt(xdim),sqrt(xdim)))
% title('reconstructed patch at the end of sampling trial')
title('b')
subplot(n_patches,3,3)
imagesc(reshape(E,sqrt(xdim),sqrt(xdim)))
% title('average reconstructed patch across sampling trial')
title('c')
colormap('gray')


figure
imagesc(squeeze(sum(S,2)))
colorbar 
title('Spike counts')
xlabel('image patch')
ylabel('neuron')

% autocorrelation for sample neurons.
acf(squeeze(S(10,:,50))',15)
title('autocorrlelation for neuron 10 on image 50')



% find PF of neurons that are mostly assembles each other
% SSIM, best is 1, worst is -1
Diff = zeros(sdim,sdim);
for i = 1:sdim
  for j = 1:sdim
    Diff(i,j) = ssim(reshape(W(:,j),8,8),reshape(W(:,i),8,8));
  end
end
[row,col] = find((abs(Diff-1)<0.5) & (abs(Diff-1)>0.00001));

for i = 1:floor(length(row)/2)
%   if((sum(sum(S(row(i),:,:))))>3) && ((sum(sum(S(col(i),:,:))))>3)
    figure
    subplot(2,1,1)
    imagesc(reshape(W(:,row(i)),8,8))
    title(strcat('neuron ', num2str(row(i))))
    colorbar
    subplot(2,1,2)
    imagesc(reshape(W(:,col(i)),8,8))
    title(strcat('neuron ', num2str(col(i))))
    colormap('gray')
    colorbar
%   end
  
end


% neuron 124 and neuron 72
[acor,lag] = xcorr(squeeze(S(124,:,50))',squeeze(S(124,:,50))',15);








% plot fano factor as 1-p, assuming each trail mean is np, variance is np(1-p)
P = sum(S,2);





% Calculating Spike Rate
% divide into 50 bins
n_patches = 50;
binsz = 20;
A = ones(1,n_sample/binsz);
Z = A*binsz;
Rate = mat2cell(S,sdim,Z,n_patches);
RateMatrix = zeros(sdim,n_sample/binsz,n_patches);
for i  = 1:length(Rate)

  rate = sum(Rate{i},2)/binsz; % dim of rate = basis * n_patches
  RateMatrix(:,i,:) = rate; 
end
Mean = squeeze(mean(RateMatrix,2));
Var = squeeze(var(RateMatrix,0,2));


%  Histogram of posterior probability distribution:
% Each neuron corresponding to each image patch
figure(2)
S_t = S(:,100:end,1);
Post =1- sum(S_t,2)/(size(S_t,2));
hist(Post(:),50)
xlabel('1-P')
ylabel('count')
title('Fano Factor distribution for all neurons across all patches')
%  set(gca,'YScale','log')



% Reconstruction error 
% histogram of average reconstruction error across all patches
figure
A = ((X - W*squeeze(S(:,end,:))).^2);
A = (sum(A,1))./(sdim);

hist(A(:),20);
title('Distribution of reconstruction error across all patches')
xlabel('average reconstruction error per pixel')
ylabel('number of image patches')

% Mean and Median as a function of samples
figure
grad = 50;
Mean_err = zeros(1,(n_sample/grad));
Median_err = Mean_err;
sample = Median_err;
for i = 1:(n_sample/grad)
A = ((X - W*squeeze(S(:,i*grad,:))).^2);
A = (sum(A,1))./(sdim);
Mean_err(i) = mean(A);
Median_err(i)=median(A);
sample(i) = i;
end
plot(sample,Mean_err,'-o')
hold on
plot(sample,Median_err,'-o')
legend('Mean Sampling Error','Median Sampling Error')

title('Reconstruction Error across Sampling Stages')
ylabel('reconstruction error per pixel')

xlabel('sampling step x50')






figure(3)
Post_m = mean(Post,3);
Post_m = squeeze(Post_m);
hist(Post_m,200)
xlabel('average posterior probability for each neuron across all image patches')

figure(4)
Post_m = mean(Post,1);
Post_m = squeeze(Post_m);
hist(Post_m,200)
xlabel('average posterior probability for each image patch across all neurons')


figure(5)
Mean = mean(S(:,100:end,:),2);
Var = var(S(:,100:end,:),0,2);
plot(Mean(:),Var(:),'+')
xlabel('Mean')
ylabel('Variance')
title('Mean and Variance of Spike Train')

figure(6)
Fano = Var(:)./Mean(:);
hist(Fano(:))
title('histogram of Fano factor distribution')




% autocorrelation for sample spike train
figure(8)


title('reconstruction error across all patches')

subplot(223)
plot(Mean(:),Var(:),'+')
xlabel('Mean')
ylabel('Variance')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('(c)')

subplot(224)
Fano = Var./Mean;
hist(Fano(:),300);
title('(d)')

print(['Data' num2str(sdim) ' basis & ' num2str(xdim) ' dimensions'],'-dpng')



% Autocorrelation function for spike train

M = Mean(:);
V = Var(:);
MM = [];
VV = [];
for i = 1:length(Mean(:))
  if(M(i) ~= 0 && V(i) ~=0)
    MM = [MM,M(i)];
    VV = [VV,V(i)];
  end
end


mat
plot(reshape(Mean(:,1,:),1,3500),reshape(Var(:,1,:),1,3500),'*')
plot(Mean, Var)

% Do first order statistics on PI

