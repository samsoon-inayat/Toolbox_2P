function data = populate_mean_raster_fitting(data,trials)

if ~exist('trials','var')
    trials = 3:10;
end

temp = load(sprintf('mean_raster_fits_file_%.2d_%.2d.mat',trials(1),trials(end)));
mrfs = temp.mean_raster_fits;

for ii = 1:6
    [rs,coeff] = getMRFS_vals(mrfs{ii});
    as = coeff(:,1);
    bs = coeff(:,2);
    cs = coeff(:,3);
    PWs = 2.36*cs./sqrt(2)*(142/50);
    data{ii}.rs = rs;
    data{ii}.coeff = coeff;
    data{ii}.pws = PWs;
    data{ii}.centers = bs * 142/50;
    data{ii}.peaks = as;
    data{ii}.sel = find(rs > 0.5 & (bs > 1 & bs < 50)' & (PWs > 5 & PWs < 120)');
    data{ii}.formula = mrfs{1}.gauss1Formula;
end