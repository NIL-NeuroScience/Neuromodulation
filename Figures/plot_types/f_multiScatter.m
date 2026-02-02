%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           f_multiScatter
% author - Brad Rauscher (created 2023)
% 
% Plots data in 'X' and 'Y' in multiple scatter plots. Plots regression
% lines for each scatter plot, and the mean regression line of all scatter
% plots (black).
% 
% INPUTS: f_multiScatter(X, Y, _)
%   X: cell array of X values (N x 1)
%   Y: cell array of Y values (N x 1)
% 
% OPTIONAL INPUTS:
%   cmp: colormap for scatter plots
%   ylabel: y-axis label
%   xlabel: x-axis label
%   title: title
%   ylim: limit for y-axis
%   xlim: limit for x-axis
%   alpha: transparency value
%   lineWidth: line width for regression lines
%   trend: plots an average trend line for all scatter plots
% 
% OUTPUTS:
%   r: correlation value(s)
%   f: figure handle
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_multiScatter(X, Y, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'cmp', [0, 0, 0]);
addParameter(p, 'ylabel', '');
addParameter(p, 'xlabel', '');
addParameter(p, 'title', '');
addParameter(p, 'ylim', []);
addParameter(p, 'xlim', []);
addParameter(p, 'alpha', 1);
addParameter(p, 'lineWidth', 1);
addParameter(p, 'trend', 1);

parse(p, varargin{:});

N = numel(X);
cIdx = round(linspace(1, size(p.Results.cmp, 1), N));

hold on;
for i = 1 : N
    scatter(X{i}, Y{i}, 'filled', ...
        MarkerFaceAlpha = p.Results.alpha, ...
        MarkerFaceColor = p.Results.cmp(cIdx(i), :));
end

% adds x limit
if ~isempty(p.Results.xlim)
    xlim(p.Results.xlim);
    current_x = xlim;
else
    current_x = xlim;
    xlim(current_x);
end

% adds y limit
if ~isempty(p.Results.ylim)
    ylim(p.Results.ylim);
else
    current_y = ylim;
    ylim(current_y);
end

% plot linear regression
ends = zeros(N, 2);
for i = 1:N
    lm = fitlm(X{i}, Y{i}); % fit linear regression to each scatter
    lm = table2array(lm.Coefficients);
    ends(i, :) = current_x * lm(2, 1) + lm(1, 1);
    plot(current_x, ends(i, :), ...
        color = [p.Results.cmp(cIdx(i),:), p.Results.alpha], ...
        lineWidth = p.Results.lineWidth);
end

% plot average trend line
if p.Results.trend
    plot(current_x, mean(ends), ...
        color = [0, 0, 0], ...
        lineWidth = 2 * p.Results.lineWidth);
end

title(p.Results.title);
xlabel(p.Results.xlabel);
ylabel(p.Results.ylabel);

end