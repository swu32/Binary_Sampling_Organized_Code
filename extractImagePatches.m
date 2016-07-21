function X = extractImagePatches( samples, winsize )

% getdata - gathers image patches from our images
%
% SYNTAX:
% getdata( samples, winsize );
%
% INPUT variables:
% samples            total number of patches to take
% winsize            patch width in pixels
%
% OUTPUT variables:
% X                  image patches in the columns
%

% How many images are there?
imagenum = 100;% 

% This will hold the patches
X=zeros(winsize^2,samples);
totalsamples = 0;

% Don't sample too close to the edges
BUFF=4;


% Load an image
Listing = dir('/Users/shuchenwu/Desktop/csc241/vanhateren_iml/*.iml');
file_list = {Listing.name};
N = 1024;

% f=-N/2:N/2-1;
% [fx fy] = meshgrid(f);
% [theta rho]=cart2pol(fx,fy);
% f_0=0.4*N;
% filtf = rho.*exp(-(rho/f_0).^4);

% Step through the images
for i=1:(imagenum)

  

% Display progress
% string = strcat('/Users/shuchenwu/Desktop/csc241/vanhateren_iml/',file_list{i})
f1 = fopen(strcat('/Users/shuchenwu/Desktop/csc241/vanhateren_iml/',file_list{i+200}),'rb','ieee-be');
w=1536;h=1024;
I=fread(f1,[w,h],'uint16');
% normalize image same as preprocessing
I = I(1:1024,1:1024);
% normalize by dividing by variance
mean_val = mean(I(:)); 
varian_val = var(I(:)); 
I = (I-mean_val)/varian_val; 
% 
% If=fftshift(fft2(I));
% 
% Iwf = filtf.*If;
% Iw = ifft2(fftshift(Iwf));
% 
% 
% D=16;
% [x y] = meshgrid(-D/2:D/2-1);
% G = exp(-0.5*((x.^2+y.^2)/(D/2)^2));
% G = G/sum(G(:));
% imv = conv2(Iw.^2,G,'same');
% imn = Iw./sqrt(imv);
% imn = I;
% imagesc(imn), axis image
% title('(d)')% constrast normalized
% colormap(gray);
% imagesc(I);


% I = double(imread(['data/' num2str(i) '.tiff']));
% I = I-mean(mean(I));
% I = I/sqrt(mean(mean(I.^2)));

% Determine how many patches to take
getsample = floor(samples/imagenum);
if i==imagenum, getsample = samples-totalsamples; end

% Extract patches at random from this image to make data vector X
for j=1:getsample
r=BUFF+ceil((size(I,1)-winsize-2*BUFF)*rand);
c=BUFF+ceil((size(I,2)-winsize-2*BUFF)*rand);
totalsamples = totalsamples + 1;
temp = reshape( I(r:r+winsize-1,c:c+winsize-1),winsize^2,1);
temp = (temp-mean(temp(:)))/(sqrt(var(temp(:))));% normalize bt standard deviation
X(:,totalsamples) =temp;
end

end

fprintf('\n');

return;

