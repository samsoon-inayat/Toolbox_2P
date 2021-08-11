function out = find_trial_by_trial_Gauss_Fit_Results(Rs,mRs,resp)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        [rsM,MFRM,centersM,PWsM] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters_trials(R.gauss_fit_on_mean,R.bin_width);
        centers(isnan(rs)) = NaN;
        jt = nanstd(centers,[],2);
        jitter_centers{rr,cc} = jt(resp{rr});
        
        gr = sum(~isnan(rs),2);
        good_rs{rr,cc} = gr(resp{rr});
        diff_centers_from_mean{rr,cc} = centers(resp{rr},:)-centersM(resp{rr})';
%         cn = 2; this_Raster = R.sp_rasters1(:,:,cn); figure(100);clf;imagesc(this_Raster);colorbar; 
%         rs(cn,:)
%         centers(cn,:)
%         n = 0;
    end
end

out.jitter_centers = jitter_centers;
out.good_rs = good_rs;
out.diff_centers_from_mean = diff_centers_from_mean;