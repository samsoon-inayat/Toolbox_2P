function selC = make_selC_struct(areCells,planeNumber,conditionsAndRasterTypes,zMI_threshold,fwidth_limits,fcenter_limits,frs_threshold,HaFD_th,HiFD_th,FR)

selC.areCells = areCells; % to see if it is identified as cell by Suite2P
selC.plane_number = planeNumber; % to select a particular plane to analyze
selC.conditionsAndRasterTypes = conditionsAndRasterTypes;
selC.zMI_threshold = zMI_threshold; % to select tuned cells
selC.fwidth_limits = fwidth_limits; % to select limits of field widths to be included
selC.fcenter_limits = fcenter_limits; % to select field center locations to be included
selC.frs_threshold = frs_threshold; % to select a threshold to be put on fitting r-square value
selC.HaFD_threshold = HaFD_th;
selC.HiFD_threshold = HiFD_th;
selC.FR_threshold = FR;
