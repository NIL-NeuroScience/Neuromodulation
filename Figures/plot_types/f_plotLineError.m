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

parse(p, varargin{:});

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
        YScale = log', ...
        XScale = log');
end

hold on;
fill([x; flipud(x)], [(y + error); flipud(y - error)], color, ...
    FaceAlpha = 0.3, ...
    EdgeColor = 'none');
plot(x, y, ...
    Color = color, ...
    LineWidth = p.Results.lineWidth);

end