function o = prop_op(prop1,prop2,th)


for rr = 1:size(prop1,1)
    for cc = 1:size(prop1,2)
        zMI_D = prop1{rr,cc}; zMI_T = prop2{rr,cc};
        diff_D_T{rr,cc} = zMI_D - zMI_T; 
        resp_diff_D_g_T{rr,cc} = diff_D_T{rr,cc} > th; resp_diff_T_g_D{rr,cc} = diff_D_T{rr,cc} < -th;
        per_resp_d_D_g_T(rr,cc) = 100*sum(resp_diff_D_g_T{rr,cc})/length((resp_diff_D_g_T{rr,cc}));
        per_resp_d_T_g_D(rr,cc) = 100*sum(resp_diff_T_g_D{rr,cc})/length((resp_diff_T_g_D{rr,cc}));
    end
end

o.resp_D_g_T = resp_diff_D_g_T;
o.resp_T_g_D = resp_diff_T_g_D;
o.resp_D_g_T_perc = per_resp_d_D_g_T;
o.resp_T_g_D_perc = per_resp_d_T_g_D;
o.diff_D_T = diff_D_T;
