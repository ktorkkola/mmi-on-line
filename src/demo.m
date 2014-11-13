data = getSatimageData('../data/satimage.data',1);
%data = getSatimageData('../data/satimage.test');
data = remmeanData(data);
% either sphere or maxnormalize!
data = sphereData(data);
% [data,maxval]=maxnormalizeData(data);
% data.spherer = eye(data.d)/maxval;

%------------------------------------------------------------
%
od=2;
clear state;
global PLOTX PLOTY; PLOTX=1; PLOTY=2;
%init = 'lda';
init = 'random';
%init = 'pca';
randn('state',3);
state.W = g_init( data, od, init);  
figure(1);whitebg('k');clf;

%------------------------------------------------------------
rand('state',1);
state.deltaW =zeros(size(state.W));
state.movie=[];
state.iteration=0;
par.momentum=0.2;
par.movie=1;
par.verbose=1;
par.plotforces = 0;
par.normalize=0; % gradient
par.batchSize = data.Ntot;
%par.batchSize = 300;

par.max_iterations = 100;
par.initialEta = 400;
par.finalEta = 50;
par.sigma=4;
state = MmiOptimizeOnline(data,par,state);

par.max_iterations = 300;
par.initialEta = 30;
par.finalEta = 1;
par.sigma = 2;
state = MmiOptimizeOnline(data,par,state);

movie(state.movie);

