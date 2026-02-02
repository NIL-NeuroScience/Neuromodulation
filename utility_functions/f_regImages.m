%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_regImages
% author - Brad Rauscher (created 2024)
% 
% Registers all images in 'imgs' to a reference allen parcellation. Uses
% individual allen parcellations in 'parcells' for each image.
% 
% INPUTS: f_regImages(imgs, refParcells, parcells, isMask)
%   imgs: cell array of images to register (N)
%   refParcells: reference allen parcellation. (H x W x 12)
%   parcells: allen parcellation for each image in 'imgs'
%   isMask: corrects for fractional values if imgs are masks
% 
% OUTPUTS:
%   reg: registered images (H x W x N)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reg = f_regImages(imgs, refParcells, parcells, isMask)

N = numel(imgs);
reg = NaN([size(refParcells.Masks, [1, 2]), N]);

for i = 1 : N
    if ~isempty(imgs{i})
        reg(:,:,i) = f_ImgReg_allen(refParcells, parcells{i}, imgs{i}, isMask);
    end
end

end