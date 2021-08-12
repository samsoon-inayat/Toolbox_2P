function out = find_trial_by_trial_Gauss_Fit_Results(Rs,mRs,resp)

for rr = 1:size(Rs,1)
    re = resp{rr};
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        [rsM,MFRM,centersM,PWsM] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters_trials(R.gauss_fit_on_mean,R.bin_width);
        centers(isnan(rs)) = NaN;
        jt = nanstd(centers,[],2);
        gr = sum(~isnan(rs),2);
        dcm = centers-centersM';
        sdc = nanmean(dcm,2);
        tre = re & sdc < 150 & sdc > 0;
        good_rs{rr,cc} = gr(tre);
        diff_centers_from_mean{rr,cc} = dcm(tre,:);
        jitter_centers{rr,cc} = jt(tre);
        dc = diff(centers,[],2);
        diff_centers{rr,cc} = dc(tre,:);
%         cn = 2; this_Raster = R.sp_rasters1(:,:,cn); figure(100);clf;imagesc(this_Raster);colorbar; 
%         rs(cn,:)
%         centers(cn,:)
%         n = 0;
    end
end

out.jitter_centers = jitter_centers;
out.good_rs = good_rs;
out.diff_centers_from_mean = diff_centers_from_mean;
out.diff_centers = diff_centers;