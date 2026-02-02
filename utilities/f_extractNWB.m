function [rfp,gfp,rfp_HD,gfp_HD,Hb,HbO,HbT,Whisking,Pupil,Accelerometer,brain_mask,vessel_mask,allen_masks,fs,mouseInfo,sessionInfo,behCam] = f_extractNWB(nwb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to extract relevant data from nwb files.
%
% INPUTS:
%    nwb: nwb file variable (first load with nwbRead from matNWB)
%
% OUTPUTS:
%    rfp: raw red fluorescence channel
%    gfp: raw green fluorescence channel
%    rfp_HD: HD-corrected red fluorescence channel
%    gfp_HD: HD-corrected green fluorescence channel
%    Hb: calculated relative change in concentration of deoxygenated HB
%    HbO: calculated relative change in concentration of oxygenated HB
%    HbT: (Hb+HbO) relative change in concentration of total HB
%    Whisking: estimated whisker movement
%    Pupil: estimated pupil diameter (relative to full eye diameter)
%    Accelerometer: mouse movement
%    brain_mask: mask of the brain exposure (Nan)
%    vesssel_mask: mask of the large pial vessels (Nan)
%    allen_masks: masks of the allen cortical regions (Nan)
%    fs: sampling rate of all measurements (Hz)
%    mouseInfo: compiled information about the mouse
%    sessionInfo: compiled information about the session
%    behCam: raw behavior camera recording (recommended not to load, large)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract optical channels 

rfp = 100*permute(nwb.acquisition.get('rfp').data.load(),[2 1 3]);
gfp = 100*permute(nwb.acquisition.get('gfp').data.load(),[2 1 3]);
rfp_HD = 100*permute(nwb.acquisition.get('rfp_HD').data.load(),[2 1 3]);
gfp_HD = 100*permute(nwb.acquisition.get('gfp_HD').data.load(),[2 1 3]);
Hb = 1e6*permute(nwb.acquisition.get('Hb').data.load(),[2 1 3]);
HbO = 1e6*permute(nwb.acquisition.get('HbO').data.load(),[2 1 3]);

rfp(isnan(rfp)) = 0;
rfp_HD(isnan(rfp_HD)) = 0;
gfp(isnan(gfp)) = 0;
gfp_HD(isnan(gfp_HD)) = 0;
Hb(isnan(Hb)) = 0;
HbO(isnan(HbO)) = 0;

HbT = Hb+HbO;

%% remove NaN pixels (zero variance error 

%% extract behavioral readouts

Whisking = nwb.processing.get('behavior').nwbdatainterface.get('Behavior').timeseries.get('Whisking').data.load();
Pupil = nwb.processing.get('behavior').nwbdatainterface.get('Behavior').timeseries.get('Pupil diameter').data.load();
Accelerometer = nwb.processing.get('behavior').nwbdatainterface.get('Behavior').timeseries.get('Accelerometer').data.load();

%% extract masks

brain_mask = permute(nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('brain_mask').image_mask.data.load(),[2 1]);
vessel_mask = permute(nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('vessel_mask').image_mask.data.load(),[2 1]);
allen_masks = double(permute(nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('allen_masks').image_mask.data.load(),[2 1 3]));
allen_masks(allen_masks==0) = NaN;

%% extract framerate

fs = nwb.acquisition.get('Hb').starting_time_rate;

%% extract mouseInfo

mouseInfo = struct;
mouseInfo.ID = nwb.general_subject.subject_id;
mouseInfo.strain = nwb.general_subject.strain;
mouseInfo.GRAB = nwb.general_subject.genotype;
mouseInfo.DOB = nwb.general_subject.date_of_birth;

surgery = strsplit(nwb.general_subject.description,': ');
mouseInfo.surgery = surgery{2};
mouseInfo.sex = nwb.general_subject.sex;

%% extract session info

sessionInfo = struct;
sessionInfo.id = nwb.identifier;

session = strsplit(nwb.identifier,'/');
sessionInfo.Mouse = session{1};
sessionInfo.Date = session{2};
sessionInfo.Run = sscanf(session{3},'Run%d');

%% extract behavior video

if nargout > 16
    behCam = permute(nwb.processing.get('behavior').nwbdatainterface.get('Behavior_Video').data.load(),[2 1 3]);
end

end