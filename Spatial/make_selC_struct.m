function selC = make_selC_struct(areCells,planeNumber,conditionNumber,rasterType,zMI_threshold,fwidth_limits,fcenter_limits,frs_threshold)

selC.areCells = areCells; % to see if it is identified as cell by Suite2P
selC.plane_number = planeNumber; % to select a particular plane to analyze
selC.conditionNumber = conditionNumber;
selC.rasterType = rasterType;
selC.zMI_threshold = zMI_threshold; % to select tuned cells
selC.fwidth_limits = fwidth_limits; % to select limits of field widths to be included
selC.fcenter_limits = fcenter_limits; % to select field center locations to be included
selC.frs_threshold = frs_threshold; % to select a threshold to be put on fitting r-square value
