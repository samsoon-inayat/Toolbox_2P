function brake_vs_no_brake_reduced

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 2 7 3 4 5 3 4 5 0 0 3 5 2 3 4 5 7 3 4 5 2 3 4 5 7];
    rasterNames = {'light22T','light22T','light22T','airOnsets22T','airOnsets22T','airOffsets22T',...
                    'airOffsets22T','airOnsets22T','airOnsets22T','airOnsets22T',...
                    'airOffsets22T','airOffsets22T','airOffsets22T','motionOnsets','motionOffsets',...
                    'light22T','light22T','airOffsets22_C','airOffsets22_C','airOffsets22_C','airOffsets22_C','airOffsets22_C'....
                    'airOnsets22_C','airOnsets22_C','airOnsets22_C','airOnsets22P','airOnsets22P','airOnsets22P','airOnsets22P','airOnsets22P'
                    };
    rasterNamesTxt = {'1-L','4-AL','6-L','2-AOn','7-AOn','2-AOff','7-AOff',...
        '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','MOn','MOff',...
        '3-A','5-A','2-Arb','3-Arb','4-Arb','5-Arb','7-Arb',...
        '3-AOnc','4-AOnc','5-AOnc','2-AOnp','3-AOnp','4-AOnp','5-AOnp','7-AOnp'};
    tic
    o = get_data(ei,selContexts,rasterNames);
    toc
    var_names = {'Lb','ArL_L','Lbs','Ab_On','Abs_On','Ab_Off','Abs_Off','Ar_On','ArL_On','Ars_On','Ar_Off','ArL_Off','Ars_Off','M_On','M_Off',...
        'Ar_L','Ars_L','Ab_Offc','Ar_Offc','ArL_Offc','Ars_Offc','Abs_Offc','Ar_Onc','ArL_Onc','Ars_Onc','Ab_Onp','Ar_Onp','ArL_Onp','Ars_Onp','Abs_Onp'};
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
magfac = mData.magfac;
n = 0;

%% after loading data go to other files for graphs and statistical testing for general graphs code go to brake_vs_no_brake_general (e.g., for making VENN diagrams)