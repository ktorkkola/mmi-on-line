function data = remmeanData(data)

%
% ver 1.5   3-14-00 
% ver 1.6  11-06-02 does not mess with the axisfactor if it exists
%

if isfield(data,'mean'),
	error('Error: Called remmeanData with existing mean!');	
end

rowX = data.X';  % change into row vectors
[Ntot,d] = size(rowX);
m=mean(rowX);
rowX = rowX-repmat( m, Ntot, 1); % remove mean

data.X = rowX';
data.mean = m;

if ~isfield(data,'axisfactor'),
    data.axisfactor = max(max(abs(rowX)));
end;
