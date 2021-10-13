function o = prop_op(prop1,prop2,th,resp)


for rr = 1:size(prop1,1)
    for cc = 1:size(prop1,2)
        zMI_D = prop1{rr,cc}; zMI_T = prop2{rr,cc};
        if exist('resp','var')
            tresp = resp{rr,cc};
            diff_D_T{rr,cc} = zMI_D(tresp) - zMI_T(tresp); 
        else
            diff_D_T{rr,cc} = zMI_D - zMI_T; 
        end
        resp_diff_D_g_T{rr,cc} = diff_D_T{rr,cc} > th; resp_diff_T_g_D{rr,cc} = -diff_D_T{rr,cc} > th;
        per_resp_d_D_g_T(rr,cc) = 100*sum(resp_diff_D_g_T{rr,cc})/length(zMI_D);
        per_resp_d_T_g_D(rr,cc) = 100*sum(resp_diff_T_g_D{rr,cc})/length(zMI_D);
    end
end

o.resp_D_g_T = resp_diff_D_g_T;
o.resp_T_g_D = resp_diff_T_g_D;
o.resp_D_g_T_perc = per_resp_d_D_g_T;
o.resp_T_g_D_perc = per_resp_d_T_g_D;
o.diff_D_T = diff_D_T;
