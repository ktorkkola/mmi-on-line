function data = getSatimageData(filename, less)

%
% ver 1.5  3-14-00 input data is in one structure
%

dat=dlmread(filename,' ');
[Ntot,d] = size(dat);

if nargin==2
   dat = dat(1:3:Ntot, :); % take only third
end
[Ntot,d] = size(dat);

Nc=6; % we know there are 6 classes !

classid=[1 2 3 4 5 7]; % class identifiers, this is external knowledge
X = [];
N = [];
J = []; prevJ = 0;
for i=1:Nc
   class = find( dat(:,d)==classid(i) );
   X = [X; dat( class, 1:d-1) ];
   N = [N length(class)];
   J = [J; prevJ+1 prevJ+length(class)];
   prevJ = prevJ + length(class);
end
d = d-1;

% data is now as row vectors

data.X = X'; % data is now in column vectors
data.J = J;
data.Ntot = Ntot;
data.N = N;
data.Nc = Nc;
data.d = d;
data.axisfactor = max(max(abs(X)));
data.labels = int2str(classid'); % classlabels as strings
