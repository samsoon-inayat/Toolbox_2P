%%
add_to_path
%% load data
clear all
clc


% % for the following recording I have to adjust the trials variable ... see
% % below
% animal_id = '171034';exp_date = '2017-11-14';recordingNumber = 1;
% ei{1} = getData(animal_id,exp_date,recordingNumber,'overwrite_behavior',0);
% ei{1}.b.air_puff_f(1) = [];
% ei{1}.b.trials = 1:20;
% ei{1}.S = fit_gauss2(ei{1});
% [ei{1}.deconv.caSigAll,ei{1}.deconv.spSigAll] = getSpikes(ei{1});

% animal_id = '170865';exp_date = '2017-08-15';recordingNumber = 1;
% ei{2} = getData(animal_id,exp_date,recordingNumber,'overwrite_behavior',0);
ee = 1;
animal_id = '172376';exp_date = '2017-11-30';recordingNumber = 1;brain_location = 'Right_Hemisphere_Location1';
ei{ee} = getData(animal_id,exp_date,recordingNumber,'brain_location',brain_location,'overwrite_behavior',1);
ei{ee}.S = fit_gauss2(ei{ee});
[ei{ee}.deconv.caSigAll,ei{ee}.deconv.spSigAll] = getSpikes(ei{ee});


