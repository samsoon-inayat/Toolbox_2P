function find_the_tonicphasic_nature

%%
si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];
Rs = o.Rs(:,si); mRs = o.mR(:,si);
% 
% si = [Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
% RsR = o.Rs(:,si);
props1 = get_props_Rs(Rs,[50,100]);


%%
an = 4; cn = 6;
resp = props1.good_FR{an,cn};
for rr = an%1:size(Rs,1)
    for cc = cn%1:size(Rs,2)
        t_Rs = Rs{rr,cc};
        spR = t_Rs.sp_rasters1;
        for cn = 1:length(resp)
            if ~resp(cn)
                continue;
            end
            raster = mean(spR(:,:,cn)); fft_raster = fft(raster);
            mag = abs(fft_raster); ang = angle(fft_raster);
            figure(1000);clf;
            subplot 131; imagesc(raster);colorbar
            subplot 132; imagesc(mag,[0 100]);colorbar
            subplot 133; imagesc(ang);colorbar
%             plot(nanmean(raster));
            pause(0.3);
        end
    end
end
