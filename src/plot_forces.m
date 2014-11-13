function plot_forces(A,dC,J,axisfactor)

% ver 1.1  01-24-01 
% ver 1.2  02-06-01 removed clf to enable subplotting
% ver 1.3  02-27-01 removed drawnow at the end 
% ver 1.4  06-04-01 provision to plot any two components of a high-dimensional set
% ver 1.5  09-27-01 added to colors and to markers to make the numbers prime
% ver 1.6  10-07-02 one-dimensional plot added

%----------------------------------------------------------------
% Check which two (or just one) dimensions to plot
global PLOTX PLOTY;
if isempty(PLOTX) & isempty(PLOTY), % if both are empty, then use dims 1 & 2
    plotx = 1;
    ploty = 2;
    clear global PLOTX PLOTY;
elseif isempty(PLOTX), % if only one is empty, assume 1-d plot
    plotx = PLOTY;
    ploty = 0;
elseif isempty(PLOTY), % if only one is empty, assume 1-d plot
    plotx = PLOTX;
    ploty = 0;
else                   % 2-d plot 
    plotx = PLOTX;
    ploty = PLOTY;
end

[od,Ntot] = size(dC);
[od,Ntot] = size(A); % TODO: assert this!
[Nc,two] = size(J);

%----------------------------------------------------------------
% Color definitions for multiple classes used both in 1-d and 2-d plots
%colors  =  [[0.5 0.2 0]', [0 0.5 0]', [0 0 0.7]', [1 0 0]', [0 0 1]', [0 1 0]', [0.5 0 0.5]', [0.7 0.2 0.2]', [1 1 0]', [0 1 1]', [1 0 1]'];
%colors  =  [[1 0 0]', [0 0 1]', [0 1 0]', [1 1 0]', [0 1 1]', [1 0 1]', [1 1 1]', [0.5 0 0]', [0 0 0.5]', [0 0.5 0]', [0.5 0.5 0]', [0 0.5 0.5]', [0.5 0 0.5]',[0.5 0.5 0.5]'];
colors  =  [[1 0 0]', [0 0 1]', [0 1 0]', [1 1 0]', [0 1 1]', [1 0 1]', [1 1 1]', ...
            [0.5 0 0]', [0 0 0.5]', [0 0.5 0]', [0.5 0.5 0]', [0 0.5 0.5]', [0.5 0 0.5]',[0.5 0.5 0.5]',...
            [0.5 0.5 1]', [1 0.5 0.5]', [0.5 1 0.5]' ]; % 17 colors (9-27-01)
[foo,cN]  = size(colors);


%----------------------------------------------------------------
% Plot of one dimension at a time
if ploty==0, % draw the one-d plot as histograms 
    
    bins = min(100,round(Ntot/20));
    edges = -axisfactor:2*axisfactor/(bins-1):axisfactor;

    histograms = zeros(Nc,bins);
    for i=1:Nc,
        histograms(i,:) = histc( A(plotx, J(i,1):J(i,2)), edges ); 
    end
    % bar(edges,histograms','stack'); % bar graph
    % plot(histograms');              % all histograms with auto coloring
    hold on;
    for i=1:Nc,                       % use a predefined coloring
      plot( edges, histograms(i,:),'color', colors(:, mod(i-1,cN)+1) ); 
    end
    hold off;
    drawnow;
else

    %----------------------------------------------------------------
    % Plot of two dimensions
    
    % Color and marker definitions for 2-d plots for multiple classes
    % markersize=8;
    markersize=5;
    Fcolors =  ['r', 'b', 'g', 'y', 'c', 'm', 'w'];
    %markers = ['o', '*', 'd', 's', 'x', '+'];
    markers = ['o', '*', 'd', 's', 'x', '+', '^', 'v', 'p', 'h','>']; % 11 markers (9-27-01)
    [foo,FcN] = size(Fcolors);
    [foo,mN]  = size(markers);
    
    % visualize forces (or just points if dC is empty)
    % clf;
    hold on;
    if ~isempty(dC),
        for i=1:Nc,
            if J(i,2)-J(i,1)>1, % quiver does not like to plot a single vector
                quiver( A(plotx,J(i,1):J(i,2)), A(ploty,J(i,1):J(i,2)), ...
                    dC(plotx,J(i,1):J(i,2)), dC(ploty,J(i,1):J(i,2)), ...
                Fcolors(mod(i-1,FcN)+1) );
            end
        end
    else
        for i=1:Nc,
            plot( A(plotx, J(i,1):J(i,2)), A(ploty, J(i,1):J(i,2)), ...
                markers(mod(i-1,mN)+1),'MarkerSize',markersize, 'MarkerEdgeColor', ...
                colors(:, mod(i-1,cN)+1)); 
        end
    end
    hold off; 
    axis( axisfactor*[ (-1.0) 1.0 (-1.0) 1.0]);
    axis square;
    drawnow;
end
