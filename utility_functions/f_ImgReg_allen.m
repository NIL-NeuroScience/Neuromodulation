%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_ImgReg_allen
% author - Brad Rauscher (created 2024)
% 
% Registers all 'img' to a reference allen parcellation. 
% 
% INPUTS: f_ImgReg_allen(refParcells, parcells, img, isMask)
%   refParcells: reference allen parcellation. (H x W x 12)
%   parcells: allen parcellation for each image in 'imgs'
%   img: 2D image to register
%   isMask: corrects for fractional values if imgs are masks
% 
% OUTPUTS:
%   reg: registered image (H x W x N)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reg, regFactor] = ...
    f_ImgReg_allen(refParcells, parcells, img, isMask)

% checks type of parcells and converts parcells array to struct
type = isstruct(parcells);

if ~type
    tmpMasks = parcells;
    tmpMasks(isnan(parcells)) = 0;
    hem = zeros(12, 1);

    for i = 1 : 12
        hem(i) = numel(bwboundaries(tmpMasks(:, :, i)));
    end
    hem = mode(hem);
    
    if hem == 2
        W = size(tmpMasks, 2);
        W = round(W / 2);

        parcells = struct;
        parcells.Masks = zeros([size(tmpMasks), 2]);

        parcells.Masks(:, 1 : W, :, 1) = tmpMasks(:, 1 : W, :);
        parcells.Masks(:, W : end, :, 2) = tmpMasks(:, W : end, :);
    else
        cBFD = bwboundaries(tmpMasks(:, :, 3));
        cBFD = mean(cBFD{1});
        cBFD = cBFD(2);
        cSSP = bwboundaries(tmpMasks(:, :, 5));
        cSSP = mean(cSSP{1});
        cSSP = cSSP(2);
        
        side = cSSP < cBFD;

        parcells = struct;
        parcells.Masks = zeros([size(tmpMasks), 2]);
        parcells.Masks(:, :, :, side + 1) = tmpMasks;
    end
end

type = 1;

% determines orientation
parcells = parcells.Masks;
L = sum(parcells(:, :, :, 1), 'all');
R = sum(parcells(:, :, :, 2), 'all');

if L && R
    hem = 2;
    side = 0;
else
    hem = 1;
    if L
        side = 0;
    else
        side = 1;
    end
end

% find points to register

if hem == 1
    Vis2 = parcells(:, :, 12, side + 1);
    SSp_ll2 = parcells(:, :, 5, side + 1);
    points.p2.point1 = round(findTop(Vis2));
    points.p2.point2 = round(findTop(SSp_ll2));
    
    Vis1 = refParcells.Masks(:, :, 12, side + 1);
    SSp_ll1 = refParcells.Masks(:, :, 5, side + 1);
    points.p1.point1 = round(findTop(Vis1));
    points.p1.point2 = round(findTop(SSp_ll1));
else
    Vis2L = parcells(:, :, 12, 1);
    Vis2R = parcells(:, :, 12, 2);
    points.p2.point1 = round(findTop(Vis2L));
    points.p2.point2 = round(findTop(Vis2R));

    Vis1L = refParcells.Masks(:, :, 12, 1);
    Vis1R = refParcells.Masks(:, :, 12, 2);
    points.p1.point1 = round(findTop(Vis1L));
    points.p1.point2 = round(findTop(Vis1R));
end


% adjust tilt

regFactor.dim1 = size(refParcells.Masks(:, :, 1, 1));

tilt(1) = atan((points.p1.point1(1) - points.p1.point2(1)) / ...
    (points.p1.point1(2) - points.p1.point2(2)));
tilt(2) = atan((points.p2.point1(1) - points.p2.point2(1)) / ...
    (points.p2.point1(2) - points.p2.point2(2)));

d_tilt = tilt(1) - tilt(2);

dim = size(parcells(:, :, 1, 1));
regFactor.tilt = d_tilt;

% rotate

cx = points.p2.point1(2);
cy = points.p2.point1(1);

angle = d_tilt * 180 / pi;
T1 = [1, 0, -cx; 0, 1, -cy; 0, 0, 1];
R = [cosd(angle), -sind(angle), 0; sind(angle), cosd(angle), 0; 0, 0, 1];
T2 = [1, 0, cx; 0, 1, cy; 0, 0, 1];

T = T2 * R * T1;

tform = affine2d(T');
img = imwarp(img, tform, OutputView = imref2d(size(img)));

% adjust scale

dist(1) = sqrt(sum((points.p1.point1 - points.p1.point2).^2));
dist(2) = sqrt(sum((points.p2.point1 - points.p2.point2).^2));

scale = dist(1) / dist(2);

xg = 1 : dim(1);
yg = 1 : dim(2);
F = griddedInterpolant({xg, yg}, img);

xq = (1 / scale : 1 / scale : dim(1))';
yq = (1 / scale : 1 / scale : dim(2))';
vq = F({xq, yq});

regFactor.scale = scale;
regFactor.dim2 = size(vq);

scaledPoint = round(points.p2.point1 * scale);

% align image sizes

reg = NaN(regFactor.dim1);

dim = size(vq);
c1 = max([1, 1; points.p1.point1 - scaledPoint + 1]);
c2 = min([regFactor.dim1; points.p1.point1 - scaledPoint + 1 + dim - 1]);

reg_c2 = c2 - points.p1.point1 + scaledPoint;
reg_c1 = scaledPoint - points.p1.point1 + c1;

reg(c1(1) : c2(1), c1(2) : c2(2)) = ...
    vq(reg_c1(1) : reg_c2(1), reg_c1(2) : reg_c2(2));

if side == 0
    reg = fliplr(reg);
end

if nargin == 4 && isMask
    reg(reg == 0) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% finds the topmost pixel in a image mask
function [point] = findTop(img)
    point = find(sum(img, 2));
    point = point(1);
    point(2) = mean(find(img(point, :)));
end

end