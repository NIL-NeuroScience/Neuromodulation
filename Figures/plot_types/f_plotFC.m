%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_plotFC
% author - Brad Rauscher (created 2024)
% 
% Plots functional connectivity matrix.
% 
% INPUTS: f_plotFC(map, diagVal, varargin)
%   map: connectivity matrix
%   diagVal: value to replace diagonal values with
% 
% OPTIONAL INPUTS:
%   cmp: colormap
%   clim: c-axis limits
%   clabel: colorbar label
%   title: title
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plotFC(map, diagVal, varargin)

% adjust diagonal values
map(diag(true(12, 1))) = diagVal;

% parse inputs
p = inputParser;
addParameter(p, 'cmp', jet);
addParameter(p, 'clim', [0,1]);
addParameter(p, 'clabel' ,'');
addParameter(p, 'title' ,'');

parse(p, varargin{:});

hold on;
imagesc(map);
axis image off;

c = colorbar;
colormap(p.Results.cmp);
clim(p.Results.clim);
c.Label.String = p.Results.clabel;

title(p.Results.title);

set(gca, FontSize = 14);

end