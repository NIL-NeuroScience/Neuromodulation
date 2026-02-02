%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_plotAllenMap
% author - Brad Rauscher (created 2024)
% 
% Plots each value in data as an allen atlas region.
% 
% INPUTS: f_plotAllenMap(data, varargin)
%   data: 12-length vector for each allen atlas region
% 
% OPTIONAL INPUTS:
%   mask: additional transparency mask
%   parcellation: allen atlas parcellation
%   cmp: colormap
%   title: title
%   cLabel: colorbar label
%   clim: c-axis limits
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_plotAllenMap(data, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'mask', []);
addParameter(p, 'parcellation', []);
addParameter(p, 'cmp', []);
addParameter(p, 'title', []);
addParameter(p, 'cLabel', []);
addParameter(p, 'clim', []);

parse(p, varargin{:});

parcellation = p.Results.parcellation;
mask = p.Results.mask;

allen_path = fullfile(f_path, 'Figures/plot_types/refAllen.mat');

if isempty(parcellation) && isempty(p.Results.mask)
    parcellation = load(allen_path);
    p.Results.mask = parcellation.refBM;
    parcellation = parcellation.refParcellation;
elseif isempty(parcellation)
    parcellation = load(allen_path);
    parcellation = parcellation.refParcellation;
end

masks = sum(parcellation.Masks, 4);

if ~isempty(mask)
    
    mask(isnan(mask)) = 0;
    masks = masks .* mask;
end

masks = logical(masks);
img = NaN(size(masks, [1, 2]));

for i = 1 : size(masks, 3)
    img(masks(:, :, i)) = data(i);
end

imAlpha = ~isnan(img);

imagesc(img, AlphaData = imAlpha);
axis image off;

% edit colormap
if ~isempty(p.Results.cmp)
    colormap(p.Results.cmp);
end

% add colorbar label
c = colorbar;
c.Label.String = p.Results.cLabel;

title(p.Results.title);
set(gca, FontSize = 14);

if ~isempty(p.Results.clim)
    clim(p.Results.clim);
end