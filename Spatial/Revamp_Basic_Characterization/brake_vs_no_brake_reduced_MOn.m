function brake_vs_no_brake_reduced

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [0 0];
    rasterNames = {'motionOnsets','motionOffsets'
                    };
    rasterNamesTxt = {'MOn','MOff'};
    tic
    o = get_data(ei,selContexts,rasterNames);
    toc
    var_names = {'M_On','M_Off'};
    for ii = 1:length(var_names)
        cmdTxt = sprintf('%s = %d;',var_names{ii},ii); eval(cmdTxt);
    end

%     Lb = 1; ArL_L = 2; Lbs = 3; Ab_On = 4; Abs_On = 5; Ab_Off = 6; Abs_Off = 7; 
%     Ar_On = 8; ArL_On = 9; Ars_On = 10; Ar_Off = 11; ArL_Off = 12; Ars_Off = 13;
%     M_On = 14; M_Off = 15; Ar_D = 16; ArL_D = 17; Ars_D = 18; Ar_T = 19; ArL_T = 20; Ars_T = 21;
%     Ar_L = 22; Ars_L = 23; Ab_Offc = 24; Ar_Offc = 25; ArL_Offc = 26; Ars_Offc = 27; Abs_Offc = 28;
%     Ar_Onc = 29; ArL_Onc = 30; Ars_Onc = 31;
%     Ab_Onp = 32; Ar_Onp = 33; ArL_Onp = 34; Ars_Onp = 35; Abs_Onp = 36;
    break
end
n = 0;

%% after loading data go to other files for graphs and statistical testing for general graphs code go to brake_vs_no_brake_general (e.g., for making VENN diagrams)