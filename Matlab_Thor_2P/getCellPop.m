function selCells = getCellPop(type)

data = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;


selCells1 = selectCells(data,mData,'Context',1);
selCells2 = selectCells(data,mData,'Context',2);
selCells3 = selectCells(data,mData,'Context',3);
selCells4 = selectCells(data,mData,'Context',4);

all = unique([selCells1 selCells2 selCells3 selCells4]);

C_1 = selectCells(data,mData,'Context',1);
    C_1_2_O = setdiff(selCells1,selCells2);
        C_1_2_O_3_O = setdiff(C_1_2_O,selCells3);
            C_1_2_O_3_O_4_O = setdiff(C_1_2_O_3_O,selCells4);
            C_1_2_O_3_O_4_P = intersect(C_1_2_O_3_O,selCells4);
        C_1_2_O_3_P = intersect(C_1_2_O,selCells3);
            C_1_2_O_3_P_4_O = setdiff(C_1_2_O_3_P,selCells4);
            C_1_2_O_3_P_4_P = intersect(C_1_2_O_3_P,selCells4);
    C_1_2_P = selectCells(data,mData,'Common',[1 2]);
        C_1_2_P_3_O = setdiff(C_1_2_P,selCells3);
            C_1_2_P_3_O_4_O = setdiff(C_1_2_P_3_O,selCells4);
            C_1_2_P_3_O_4_P = intersect(C_1_2_P_3_O,selCells4);
        C_1_2_P_3_P = intersect(C_1_2_P,selCells3);
            C_1_2_P_3_P_4_O = setdiff(C_1_2_P_3_P,selCells4);
            C_1_2_P_3_P_4_P = intersect(C_1_2_P_3_P,selCells4);
    C_2 = setdiff(selCells2,C_1);
        C_2_3_O = setdiff(C_2,selCells3);
            C_2_3_O_4_O = setdiff(C_2_3_O,selCells3);
            C_2_3_O_4_P = intersect(C_2_3_O,selCells3);
        C_2_3_P = intersect(C_2,selCells3);
            C_2_3_P_4_O = setdiff(C_2_3_P,selCells3);
            C_2_3_P_4_P = intersect(C_2_3_P,selCells3);
        C_3 = setdiff(selCells3,selCells2);
            C_3_4_O = setdiff(C_3,selCells4);
            C_3_4_P = intersect(C_3,selCells4);
            C_3_1_P = intersect(C_3,selCells1);
            C_4 = setdiff(selCells4,selCells3);
        
    C_2a = selectCells(data,mData,'Context',2);
        C_2a_3_O = setdiff(C_2a,selCells3);
            C_2a_3_O_4_O = setdiff(C_2a_3_O,selCells3);
            C_2a_3_O_4_P = intersect(C_2a_3_O,selCells3);
        C_2a_3_P = intersect(C_2a,selCells3);
            C_2a_3_P_4_O = setdiff(C_2a_3_P,selCells3);
            C_2a_3_P_4_P = intersect(C_2a_3_P,selCells3);
        
        C_3a = selectCells(data,mData,'Context',3);
            C_3a_4_O = setdiff(C_3a,selCells4);
            C_3a_4_P = intersect(C_3a,selCells4);

            C_4a = selectCells(data,mData,'Context',4);
% 
% if strcmp(type,'cells of context 1 disrupted in context 2')
%     selCells = C_1_2_O;
%     return;
% end
% 
% if strcmp(type,'cells of context 1 that did not get disrupted in context 2')
%     selCells = C_1_2_P;
%     return;
% end
% 

cmdTxt = sprintf('selCells = %s;',type);
eval(cmdTxt);
n = 0;