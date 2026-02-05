%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           f_plotLineError
% author - Brad Rauscher (created 2024)
% 
% Plots line with errorbars.
% 
% INPUTS: f_plotLineError(x, y, error, varargin)
%   x: x values
%   y: y values
%   error: errorbar values to +/- to y
%   xlim: x-axis limits
%   xlim: y-axis limits
%   xlabel: x-axis label
%   ylabel: y-axis label
%   title: title
%   legend: cell array of legend entries
% 
% OPTIONAL INPUTS:
%   color: line color
%   lineWidth: line width
%   log: sets axis scale to log
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plotLineError(x, y, error, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'color', []);
addParameter(p, 'lineWidth', 2);
addParameter(p, 'log', 0);
addParameter(p, 'xlim', []);
addParameter(p, 'ylim', []);
addParameter(p, 'xlabel', []);
addParameter(p, 'ylabel', []);
addParameter(p, 'title', []);
addParameter(p, 'legend', []);

parse(p, varargin{:});

x = x(:);
y = y(:);
error = error(:);

plotLegend = p.Results.legend;
if iscell(plotLegend)
    plotLegend = [plotLegend; repmat({''}, 1, numel(plotLegend))];
    plotLegend = plotLegend(:);
end

if isempty(p.Results.color)
    color = get(groot, 'defaultAxesColorOrder');
    color = color(1, :);
else
    color = p.Results.color;
end

% set axis to log and remove 0 x values
if p.Results.log
    zeroIdx = x == 0;
    x(zeroIdx) = [];
    y(zeroIdx) = [];
    error(zeroIdx) = [];
    set(gca, ...
        YScale = 'log', ...
        XScale = 'log');
end

hold on;
fill([x; flipud(x)], [(y + error); flipud(y - error)], color, ...
    FaceAlpha = 0.3, ...
    EdgeColor = 'none');
plot(x, y, ...
    Color = color, ...
    LineWidth = p.Results.lineWidth);

if ~isempty(p.Results.xlim)
    xlim(p.Results.xlim);
end

if ~isempty(p.Results.ylim)
    ylim(p.Results.ylim);
end

xlabel(p.Results.xlabel);
ylabel(p.Results.ylabel);
title(p.Results.title);
if ~isempty(plotLegend)
    legend(plotLegend);
end

end