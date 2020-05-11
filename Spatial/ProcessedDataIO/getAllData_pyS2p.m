function ei = getAllData_pyS2p(animal_FS,animal_id,exp_date,recordingNumber,brain_location,ow)
if ~exist('ow','var')
    ow = 0;
end
% animal_id = '171179';exp_date = '2018-02-12';recordingNumber = 2;brain_location = 'Right_Hemisphere_Location1';
ei = getData_pyS2p(animal_FS,animal_id,exp_date,recordingNumber,'brain_location',brain_location,'overwrite_behavior',ow);
% ei.S = fit_gauss2(ei);
[ei.deconv.caSigAll,ei.deconv.spSigAll] = getSpikes(ei);

