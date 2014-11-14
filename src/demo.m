
clear; close all;
data = getSatimageData('../data/satimage.data',1);
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
randn('state',1);
state.W = g_init( data, od, init);  
figure(1);whitebg('k');clf;

%------------------------------------------------------------
%rand('state',1);
state.deltaW =zeros(size(state.W));
state.movie=[];
state.iteration=0;
par.momentum=0.2;   % on-line learning momentum
par.movie=1;        % generate a movie or no?
par.verbose=1;      % text to console
par.plotforces = 0; % MMI evaluated only if this is 1
par.normalize=0;    % normalize gradient or not
% batch_size: how many random pairs of points taken. 
% Keep as large as possible
par.batchSize = 10*data.Ntot; 
%par.batchSize = 500; % smaller batch requires more iterations!

par.max_iterations = 200;
par.initialEta = 100; % learning rate 
par.finalEta = 100;
par.sigma=3.0;   % Width of the Gaussian kernel. Important! Depends on data.
state = MmiOptimizeOnline(data,par,state);

par.max_iterations = 200;
par.initialEta = 50;
par.finalEta = 10;
par.sigma=2.0;
state = MmiOptimizeOnline(data,par,state);

movie(state.movie);

