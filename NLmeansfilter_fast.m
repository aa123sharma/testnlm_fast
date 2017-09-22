% A few modifications done in the implementation by the original author Jose Vicente Manjon Herrera
... 1) Enabled pca based similarity weight computation
... 2) Moved all the loop computations inside a mex file called NLM_fast*
... 3) Similarity is now computed using patches per pixel matrix V, which is obtained from another mex file called image2_vectors_double*
...    (This is from the very efficient implementation of NLM by Dirk-Jan Kroon. Code -  https://www.mathworks.com/matlabcentral/fileexchange/27395-fast-non-local-means-1d--2d-color-and-3d
... 4) For pca based similarity, added pca function here
... 5) image2_vectors_double already returns gaussian smoothed patches, so removed 'kernel' smoothing 
... 6) Although not done here, the function NLM_fast can only be used to compute the Weight matrix (defining weights per pixel in a N_totalxN_total sparse matrix)
... 7) Additional support for RGB images
    
... 8)The implementation by Dirk-Jan Kroon is still the best, but I needed a simple and efficient code for also generating the Weight matrix

function [output, W, W_norm] = NLmeansfilter_fast(input,t,f,h,pca_en,pca_ne)
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  input: image to be filtered
 %  t: radio of search window
 %  f: radio of similarity window
 %  h: degree of filtering
 %
 %  Author: Jose Vicente Manjon Herrera & Antoni Buades
 %  Date: 09-03-2006
 %
 %  Implementation of the Non local filter proposed for A. Buades, B. Coll and J.M. Morel in
 %  "A non-local algorithm for image denoising"
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Size of the image
 [m,n,ndims]=size(input);
 
 % Replicate the boundaries of the input image
 input2 = padarray(input,[f f],'symmetric');
 
 ... Default pca_en = 0
 if(nargin<5)
     pca_en = 0;
 end
 if(pca_en == 1)
     fprintf('Using PCA coeffcients for patch similarity!\n');
     if(nargin<6)
         error('No. of Principle Components missing!');
     end
     % Get the local patches of every pixel-coordinate in the block
     V = image2vectors_double(double(input2), double(f), 1); 
     % Do PCA on the block (Since, we are using the PCA coefficients, no need to project back)
     ... Restrict PCA coefficients to pca_ne only. V_coeff should be a (2*f+1)^2 x pca_ne matrix
     % Do PCA on the block
     [~, Evectors, x_mean]=PCA(V(:,1:end),pca_ne);
     % Project the block to the reduced PCA feature space. V should be pca_ne x (m*n) matrix
     V = Evectors'*(V-repmat(x_mean,1,size(V,2)));
     ... Convert to form (m*n)xpca_ne format
     V = V';
 else
     % Get the local patches of every pixel-coordinate in the block
     V = image2vectors_double(double(input2), double(f), 1); 
     ... Convert to form (m*n)x(2*f+1)^2*(ndims) format
     V = V';
 end

 
% Generate NLM output 
output = NLM_fast(V, input, input2, h, f, t);


function [Evalues, Evectors, x_mean]=PCA(x,ne)
% PCA using Single Value Decomposition
% Obtaining mean vector, eigenvectors and eigenvalues
%
% [Evalues, Evectors, x_mean]=PCA(x,ne);
%
% inputs,
%   X : M x N matrix with M the trainingvector length and N the number
%              of training data sets
%   ne : Max number of eigenvalues
% outputs,
%   Evalues : The eigen values of the data
%   Evector : The eigen vectors of the data
%   x_mean : The mean training vector
%
s=size(x,2);

% Calculate the mean 
x_mean=sum(x,2)/s;

% Substract the mean
x2=(x-repmat(x_mean,1,s))/ sqrt(s-1);

% Do the SVD 
[U2,S2] = svds(x2,ne,'L',struct('tol',1e-4));
Evalues=diag(S2).^2;
Evectors=U2;
