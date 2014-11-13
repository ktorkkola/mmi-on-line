function val = Gaussian(vector, covar);

% n-dimensional Gaussian with covariance 'covar'
% which is either a scalar, a row vector, or a matrix
% 
% vector - column vector
% val    - row vector result
% 
% 1.0 99-10-27 original 
% 1.1 99-11-01 revised kt.
% 1.2 00-09-11 bug correcred
% 2.0 01-03-12 all covariance types
% 2.1 01-03-25 all covariance types (really)
% 2.2 01-10-22 vector can be a matrix (in the spherical case)
% 2.3 02-11-13 corrected scalar case for 1-d vector

[dim,N] = size(vector);
[m,n]=size(covar);

if m*n==1 % scalar, spherical covariance
   
  w = -1/(2*covar);
  constant = (-w/pi)^(dim/2) ;
  val = exp( sum(vector.*vector,1) * w ) * constant;
  
elseif m==1 | n==1 % vector, diagonal covariance
      
  constant = 1 / ( sqrt(prod(covar)) * (2*pi)^(dim/2) );  
  val = exp( -0.5*sum(vector.*vector./covar') ) * constant;    
      
else % full covariance
   
  constant = 1 / ( abs(det(covar)) * (2*pi)^(dim/2) );  
  val = exp( -0.5*vector'*inv(covar)*vector ) * constant;    
  %disp('Error: not implemented yet');    
end

