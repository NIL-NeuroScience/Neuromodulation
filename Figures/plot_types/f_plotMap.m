%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_plotMap
% author - Brad Rauscher (created 2024)
% 
% Plot image. 
% 
% INPUTS: f_plotMap(map, varargin)
%   map: image to plot
% 
% OPTIONAL INPUTS:
%   cmp: colormap
%   clim: c-axis limits
%   clabel: colorbar label
%   title: title
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plotMap(map, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'cmp', jet);
addParameter(p, 'clim', prctile(map(:), [1, 99]));
addParameter(p, 'clabel', '');
addParameter(p, 'title', '');

parse(p, varargin{:});

imAlpha = ~isnan(map);

hold on;

imagesc(map, AlphaData = imAlpha);
axis image off;

c = colorbar;
colormap(p.Results.cmp);
clim(p.Results.clim);
c.Label.String = p.Results.clabel;

title(p.Results.title);

set(gca, ...
    FontSize = 14, ...
    YDir = 'reverse');

end