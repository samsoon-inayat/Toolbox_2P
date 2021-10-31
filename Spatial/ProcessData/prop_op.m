function o = prop_op(prop1i,prop2i,thi,resp)

prop1 = prop1i.zMI;
prop2 = prop2i.zMI;

for rr = 1:size(prop1,1)
    for cc = 1:size(prop1,2)
        zMI_D = prop1{rr,cc}; zMI_T = prop2{rr,cc};
        if exist('resp','var')
            tresp = resp{rr,cc};
            diff_D_T{rr,cc} = zMI_D(tresp) - zMI_T(tresp); 
        else
            diff_D_T{rr,cc} = zMI_D - zMI_T; 
        end
        if thi < 0
            th = nanmean(diff_D_T{rr,cc}) + 1*nanstd(diff_D_T{rr,cc});
        else
            th = thi;
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
% o.good_FR_D_g_T = cell_list_op(resp_diff_D_g_T,prop1i.good_FR,'and');
% o.good_FR_T_g_D = cell_list_op(resp_diff_T_g_D,prop2i.good_FR,'and');
