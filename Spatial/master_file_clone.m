%% SET ALL DEFAULT OPTIONS HERE

% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter

% check out the README file for detailed instructions (and extra options)
% addpath('D:\CODE\MariusBox\runSuite2P') % add the path to your make_db file

% overwrite any of these default options in your make_db file for individual experiments
% make_db; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE

% ops0.toolbox_path = 'G:\Sam\OneDrive - University of Lethbridge\Google Drive\MatlabCode\Suite2P';
% if exist(ops0.toolbox_path, 'dir')
% 	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
% else
% 	error('toolbox_path does not exist, please change toolbox_path');
% end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION
% if strcmp(cName,'Sam-WS')
%     ops0.useGPU                 = 1; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
% else
%     ops0.useGPU                 = 0; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
% end
ops0.useGPU                 = 0;
ops0.fig                    = 1; % turn off figure generation with 0
% ops0.diameter               = 20; % most important parameter. Set here, or individually per experiment in make_db file

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            = makeName('Tif_Data',mainFolder); % Suite2P assumes a folder structure, check out README file
ops0.temp_tiff              = makeName('temp1.tif',mainFolder); % copies each remote tiff locally first, into this file
ops0.RegFileRoot            = makeName('Processed_Data',mainFolder);  % location for binary file
ops0.DeleteBin              = 0; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = makeName('Processed_Data',pmainFolder); % a folder structure is created inside
ops0.RegFileTiffLocation    = '';%makeName('Processed_Data',mainFolder); %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)

% registration options
ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.dobidi                 = 0; % infer and apply bidirectional phase offset
ops0.nonrigid               = 0;

% cell detection options
ops0.getROIs                = 1; % whether to run optimization for finding ROIs
ops0.ShowCellMap            = 1; % during optimization, show a figure of the clusters
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 1000; % how many SVD components for cell clustering
ops0.signalExtraction       = 'raw'; % how to extract ROI and neuropil signals: 'raw', 'regression'

% spike deconvolution options

ops0.sensorTau              = 2; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = Inf; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)

% red channel options
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

%% RUN THE PIPELINE HERE
db0 = db;

for iii = 1:length(dbii) %[3:length(db) 1:2]
    iexp = dbii(iii);
    ops0.imageRate              = ei{iii}.frameRate;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
    ops0.NavgFramesSVD          = floor(ei{iii}.totalFrames/4); % how many (binned) timepoints to do the SVD based on
    ops0.NimgFirstRegistration  = floor(ei{iii}.totalFrames/40); % number of images to include in the first registration pass 
    ops0.nimgbegend             = floor(ei{iii}.totalFrames/80); % frames to average at beginning and end of blocks
    run_pipeline(db(iexp), ops0);
    
    % deconvolved data into (dat.)cl.dcell, and neuropil subtraction coef
%     add_deconvolution(ops0, db0(iexp));
%     do_deconvolution(ops0,db0(iexp));
%     % add red channel information (if it exists)
%     if isfield(db0,'expred') && ~isempty(db0(iexp).expred)
%         ops0.nchannels_red = db0(iexp).nchannels_red;
%         
%         run_REDaddon(iexp, db0, ops0) ; % create redcell array
%         
%         DetectRedCells; % fills dat.cl.redcell and dat.cl.notred
%     end
    
end
%% STRUCTURE OF RESULTS FILE
% cell traces are in dat.Fcell
% neuropil traces are in dat.FcellNeu
% manual, GUI overwritten "iscell" labels are in dat.cl.iscell
%  
% stat(icell) contains all other information:
% iscell: automated label, based on anatomy
% neuropilCoefficient: neuropil subtraction coefficient, based on maximizing the skewness of the corrected trace (ICA)
% st: are the deconvolved spike times (in frames)
% c:  are the deconvolved amplitudes
% kernel: is the estimated kernel
