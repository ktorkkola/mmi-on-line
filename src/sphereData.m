function data = sphereData(data)

%
% ver 1.5   3-14-00 
% ver 1.51  2-06-01 More verbose varning.
%

% check that the mean has been removed, if not,
% remove the mean.

if ~isfield(data,'mean')
	data=remmeanData(data);
end

[data.X, newdimension, data.spherer] = sphere(data.X);
if (newdimension<data.d)
	disp(sprintf('SphereData: Low rank warning: Reduced dimension from %d to %d\n',data.d,newdimension));
end
data.d = newdimension;

% appropriate display factor after sphering
data.axisfactor = 3.0; 
