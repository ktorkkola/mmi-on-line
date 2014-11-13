function [X,d,spherer] = sphere(X)

% takes a column vector matrix X, spheres the vectors and 
% returns a column vector matrix
%
% ver 1.5  3-14-00 
%
epsilon = 1e-8;

% Change X into row vector matrix
rowX = X';
[Ntot,d] = size(rowX);

% Assume that the mean has already been removed 
% but check it anyway!
% m=mean(rowX);
% if sum(m)>epsilon
% 	rowX = rowX-repmat( m, Ntot, 1); % remove mean
%	disp('Warning: Called sphere without having removed the mean!');	
% end

covarianceMatrix = cov(rowX);
maxLastEig = rank(covarianceMatrix);		% rank of covariance
options.disp=0;
[E,D] = eigs(covarianceMatrix, maxLastEig, 'LM',options);
spherer = (inv (sqrt (D)) * E');
X =  spherer * rowX';			% sphere data
% data is now in column vectors

d = maxLastEig;
