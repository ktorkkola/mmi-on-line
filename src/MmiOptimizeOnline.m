function state = MmiOptimizeOnline(data,par,state)

% ver 1.0  08-23-01 kt 
% ver 1.1  10-22-01 Added momentum and movie.
% ver 1.2  10-22-01 Experimenting with vectorizing.
% ver 1.3  10-23-01 if batchSize <= 0, then uses all the samples
% ver 1.4  11-13-02 uses now par.verbose

%------------------------------------------------------------
if ~isfield(state,'movie')  movie=0; else movie=1; end
iteration = state.iteration;
max_iterations=par.max_iterations;
batchSize = par.batchSize;
initialEta = par.initialEta;
finalEta = par.finalEta;
verbose=par.verbose;
[od,d] = size(state.W);
I=0.0;


if ((iteration==0) & (movie==1))
   clear state.movie;
end

if isfield(par,'momentum')
    momentum = par.momentum;
    if isfield(state,'deltaW')
        deltaW = state.deltaW;
    else
        deltaW = zeros(size(state.W));
    end
else
    momentum = 0;
end


%------------------------------------------------------------
label = zeros(1,data.Ntot);
for class=1:size(data.J,1),
    for i=data.J(class,1):data.J(class,2),
        label(i)=class;
    end
end

%------------------------------------------------------------
updateW = zeros(size(state.W));
y = state.W * data.X;

%---------------------------------------------------------
% continue previous iterations if iteration is not zeroed.
bgnit=iteration+1; endit=iteration+max_iterations;
for iteration=bgnit:endit
    
    % linearly decreasing plasticity
    eta = initialEta + ((iteration-bgnit)/(endit-bgnit))*(finalEta-initialEta); 
    
    % generate indices for random samples
    if batchSize >= data.Ntot,
        % use all samples at once if possible
        i1=1:data.Ntot;
        if batchSize > data.Ntot,
            i1(data.Ntot:batchSize) = floor(rand(1,batchSize-data.Ntot+1) * data.Ntot) + 1;
        end
    else
        % take batchSize random pairs 
        i1 = floor(rand(1,batchSize) * data.Ntot) + 1;
    end    
    i2 = floor(rand(1,batchSize) * data.Ntot) + 1;
    
    ydiff = y(:,i2) - y(:,i1);
    gdiff = Gaussian(ydiff,2*par.sigma^2);
    sameClass = ( label(i1)==label(i2) );
    sameClass = (sameClass-1)*(1/8);
        
    ydiff = ydiff .* repmat( (sameClass.*gdiff),od,1);
    updateW = ((data.X(:,i1)-data.X(:,i2)) * ydiff')';

    if par.normalize,
        updateW = updateW / sqrt(sum(sum(updateW.*updateW))); % normalized gradient
    else
        updateW = updateW / (batchSize); % not needed if normalized
    end
    
       
    % momentum (add a fraction of the previous change)
    if momentum==0,
      % update and orthonormalize W
      state.W = -orthonormalize(state.W + eta*updateW);
    else
      deltaW = eta * updateW + momentum * deltaW;
      newW = -orthonormalize(state.W + deltaW);  % orthogonalize
      deltaW = newW-state.W;
      state.W = newW;
    end
    updateW = zeros(size(state.W));
    
    % project data
    y = state.W*data.X;
    
    % draw
    if ((movie==1) | (verbose>=3)),
      if isfield(par,'forcesfig'),
         figure(par.forcesfig);
      end
      clf;
      if par.plotforces==1,
         [I,dC,ma,av] = pfMex(y, par.sigma, data.J);
         plot_forces(y,dC,data.J,data.axisfactor);
      else
         plot_forces(y,[],data.J,data.axisfactor);
      end
    end
    if movie==1,
       if iteration==1,
           state.movie = getframe;
       else
           state.movie(iteration) = getframe;
       end
    end

    % display
    if (verbose>=1),
       disp(sprintf('i=%3d  MI=%f  eta=%f  sigma=%f',iteration,I,eta,par.sigma));
    end
end

state.finalEta = finalEta;
state.iteration = iteration;
state.deltaW = deltaW;


