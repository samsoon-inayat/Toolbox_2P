function R1 = combine_rasters(R)

R1 = R{1};
R1.sp_rasters1 = cat(1,R1.sp_rasters1,R{2}.sp_rasters1);
for ii = 3:length(R)
    tR = R{ii};
    R1.sp_rasters1 = cat(1,R1.sp_rasters1,tR.sp_rasters1);
end