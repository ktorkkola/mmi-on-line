function W = pca(X, od)
%
% W = pca(X,od)
%
% X  - input data as column vectors
% od - desired output space dimension
% W  - returns the basis vectors as rows
%
%
% Kari Torkkola, 000314
%                010429

% covariance of data
covariance=cov(X');

% compute the desired number of eigenvectors corresponding to largest 
% eigenvalues of the covariance matrix
d = size(covariance,1);
if d<od howmany = d; else howmany = od; end;

options.disp=0;
[eigenvectors,eigenvalues] = eigs(covariance,howmany,'lm',options);

% returned eigenvectors have length=1
%W = eigenvectors';

% scale them with eigenvalues to get variances in these directions
%W = (eigenvectors*eigenvalues)';  

% scale them with inverse eigenvalues 
%W = (eigenvectors* diag(1 ./ diag(eigenvalues)))';

% now rows of W are the directions

if howmany<od
   W = [eigenvectors rand(d, od-howmany)-0.5 ];
   if od<=d % if requested less than original dim we can orthogonalize
      W = orth(W);
   end
   W = W';
else
   W = eigenvectors';
end
