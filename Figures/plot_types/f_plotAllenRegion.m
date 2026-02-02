%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          f_plotAllenRegion
% author - Brad Rauscher (created 2024)
% 
% Plots each value in data as an allen atlas region.
% 
% INPUTS: f_plotAllenRegion(region, side, varargin)
%   region: region to plot (1 through 12)
%   side: hemisphere to plot (1 - left, 2 - right)
% 
% OPTIONAL INPUTS:
%   mask: additional transparency mask
%   parcellation: allen atlas parcellation
%   color: outline color
%   lineWidth: outline width
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plotAllenRegion(region, side, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'mask', []);
addParameter(p, 'parcellation',[ ]);
addParameter(p, 'color', [0, 0, 0]);
addParameter(p, 'lineWidth', 1);

parse(p, varargin{:});

parcellation = p.Results.parcellation;
mask = p.Results.mask;

allen_path = fullfile(f_path, 'Figures/plot_types/refAllen.mat');

if isempty(parcellation) && isempty(mask)
    parcellation = load(allen_path);
    mask = parcellation.refBM;
    parcellation = parcellation.refParcellation;
elseif isempty(parcellation)
    parcellation = load(allen_path);
    parcellation = parcellation.refParcellation;
end

if isempty(side)
    masks = sum(parcellation.Masks(:, :, region, :), 4);
else
    masks = parcellation.Masks(:, :, region, side);
end

if ~isempty(mask)
    mask(isnan(mask)) = 0;
    masks = masks .* mask;
end

bound = bwboundaries(masks);
N = numel(bound);

hold on;
for i = 1 : N
    plot(bound{i}(:, 2), bound{i}(:,1), ...
        Color = p.Results.color, ...
        LineWidth = p.Results.lineWidth);
end

set(gca, 'YDir', 'reverse');