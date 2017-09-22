% Please make sure to build the two mex files first before running this demo
... mex NLM_fast.c 
... mex image2vectors_double.c

%% Data 2-D
I1 = im2double(rgb2gray(imread('im6.png')));
In1 = zeros(size(I1));
sigman = 25/255;
for ch=1:size(I1,3)
    In1(:,:,ch) = imnoise(I1(:,:,ch), 'gaussian', 0, (sigman)^2);
end
% Set parameters
p_rad = 3;
n_rad = 9;
h     = 0.05; 
%  With NLmeansfilter
fprintf('Checking gray image with original implementation\n');
tic;
In1_out1_1 = NLmeansfilter(In1, n_rad, p_rad, h);
toc;
%  With NLmeansfilter_fast
fprintf('Checking gray image with faster implementation\n');
tic;
In1_out2_1 = NLmeansfilter_fast(In1, n_rad, p_rad, h, 0, 0);
toc;
figure(1), 
subplot(2,3,1), imshow(In1), title(sprintf('[%f]', psnr(In1,I1)));
subplot(2,3,2), imshow(In1_out1_1), title(sprintf('[%f]', psnr(In1_out1_1,I1)));
subplot(2,3,3), imshow(In1_out2_1), title(sprintf('[%f]', psnr(In1_out2_1,I1)));

%% Data 3-D
I1 = im2double((imread('im6.png')));
In1 = zeros(size(I1));
sigman = 55/255;
for ch=1:size(I1,3)
    In1(:,:,ch) = imnoise(I1(:,:,ch), 'gaussian', 0, (sigman)^2);
end
% Set parameters
p_rad = 3;
n_rad = 9;
h     = 0.2; 
%  With NLmeansfilter (This function cannot work on RGB images, hence the denoising has to be carried out per channel, which is no doubt incorrect.
... But doing this for demonstration purposes. 
fprintf('Checking RGB image with original implementation\n');
tic;
In1_out1_2 = zeros(size(In1));
for ch=1:3
    In1_out1_2(:,:,ch) = NLmeansfilter(In1(:,:,ch), n_rad, p_rad, h);
end
toc;
%  With NLmeansfilter_fast
fprintf('Checking RGB image with faster implementation\n');
tic;
In1_out2_2 = NLmeansfilter_fast(In1, n_rad, p_rad, h, 0, 0);
toc;
figure(1), 
subplot(2,3,4), imshow(In1), title(sprintf('[%f]', psnr(In1,I1)));
subplot(2,3,5), imshow(In1_out1_2), title(sprintf('[%f]', psnr(In1_out1_2,I1)));
subplot(2,3,6), imshow(In1_out2_2), title(sprintf('[%f]', psnr(In1_out2_2,I1)));
