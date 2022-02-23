function [oR1,xss,xse] = combine_rasters(R)

for an = 1:size(R,1)
    tRs = R(an,:);
    R1 = tRs{1};
    if an == 1
        xss(1) = 1; xse(1) = size(R1.sp_rasters1,2);
    end
    R1.sp_rasters1 = cat(2,R1.sp_rasters1,tRs{2}.sp_rasters1);
    if an == 1
        xss(2) = xse(1) + 1; xse(2) = size(R1.sp_rasters1,2);
    end
    for ii = 3:length(tRs)
        tR = tRs{ii};
        R1.sp_rasters1 = cat(2,R1.sp_rasters1,tR.sp_rasters1);
        if an == 1
            xss(ii) = xse(ii-1) + 1; xse(ii) = size(R1.sp_rasters1,2);
        end
    end
    oR1{an,1} = R1;
end