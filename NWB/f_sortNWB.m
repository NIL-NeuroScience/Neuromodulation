%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_sortNWB
% author - Brad Rauscher (created 2024)
% 
% Organizes .nwb files in 'path_NWB'
% 
% INPUTS: f_sortNWB(path_NWB)
%   path_NWB: path to nwb file directory
% 
% OUTPUTS:
%   nwb_list: list and information for each .nwb file in path_NWB
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nwb_list = f_sortNWB(path_NWB)

files = dir(path_NWB);
files(1 : 2) = [];

nwb_list = struct( ...
    Mouse = [], ...
    Date = [], ...
    GRAB = [], ...
    Run = [], ...
    Path = []);

N = numel(files);

for i = 1 : N
    name = fullfile(files(i).folder, files(i).name);
    nwb_list(i).Path = name;
    nwb_list(i).GRAB = h5read(name, '/general/subject/genotype');
    name = h5read(name, '/identifier');
    name = strsplit(name, '/');

    nwb_list(i).Mouse = name{1};
    nwb_list(i).Date = name{2};
    nwb_list(i).Run = str2double(name{3}(4 : 5));
end

end