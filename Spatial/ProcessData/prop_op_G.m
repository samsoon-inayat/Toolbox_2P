function o = prop_op_G(prop1i,prop2i,propN)

cmdTxt = sprintf('prop1 = prop1i.%s;',propN);eval(cmdTxt);
cmdTxt = sprintf('prop2 = prop2i.%s;',propN);eval(cmdTxt);

for rr = 1:size(prop1,1)
    for cc = 1:size(prop1,2)
        zMI_D = prop1{rr,cc}; zMI_T = prop2{rr,cc};
        if isempty(zMI_D)
            zMI_D = NaN(size(zMI_T));
        end
        diff_D_T{rr,cc} = zMI_D - zMI_T; 
        diff_T_D{rr,cc} = zMI_T - zMI_D; 
    end
end

o.diff_D_T = diff_D_T;
o.diff_T_D = diff_T_D;
