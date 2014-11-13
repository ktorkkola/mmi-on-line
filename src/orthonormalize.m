function oW = orthonormalize(iW)

%
% orthonormalizes rows of W using qr
%
% ver 1.0  4-26-00 
% 

[r,c] = size(iW);

% orthogonalize (qr orthogonalizes columns)
[oW,R] = qr(iW',0);
oW = -oW';  % sign change for convenience 

% normalize
for i=1:r
   oW(i,:) = oW(i,:) / sqrt( oW(i,:)*oW(i,:)' );
end
