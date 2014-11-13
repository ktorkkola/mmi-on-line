function [W,W0,theta,td] = g_init( data, od, init_type, initW)
%
% Initialize rotations for linear feature transformation
%
% ver 1.1  12-21-99  to accommodate the correct number of parameters od*(d-od)
% ver 1.5  3-13-00 input data is in one structure
% ver 1.6  3-18-00 added pca initialization
% ver 1.7  3-21-00 pca also with sphered data
% ver 1.8  9-20-01 quick rand init
% ver 1.9  10-9-02 initial W can be given
%
% od      - output space dimension
% data.d  - input space dimension
% td      - how many free rotations there are 
%
% W0      - initial column basis vectors...
% theta   - ...rotated by angles theta...
% W       - ...to give row basis vectors W

if nargout==1, quick=1; else quick=0; end

td = od*(data.d-od);

% initial W has been given
if nargin>3 & size(initW,1)>0,
    W = initW;
    W0 = W';
    theta = zeros(1,td);            % zero rotation angles
    return;
end

switch init_type

case 'random'
    if quick,
        W = randn(data.d,od);    
        W = orth(W)';
    else
        theta = pi*(rand(1,td)-0.5);    % initial rotation angles - the parameters
        W0=eye(data.d); W0=W0(:,1:od);  % column basis vectors - start with unit vectors
        W = g_rotate(W0,theta);         % rotated row basis vectors
    end
                                
case 'lda'
   W = lda(data.X, data.J, od);    % initialize with LDA row vectors
   W0 = W';
   theta = zeros(1,td);            % zero rotation angles
 
case 'pca'
   % if data is sphered, then it is pca'd. Take the first od components!
   if isfield(data,'spherer')
      W=eye(data.d); W=W(1:od,:);  % row basis vectors
   else
      W = pca(data.X, od);         % initialize with PCA row vectors
   end
   W0 = W';
   theta = zeros(1,td);            % zero rotation angles
 
otherwise
   error('Unknown initialization!');

end
