function zMI = get_average_zMI(R)

if isstruct(R)
    zMI = nanmean(R.info_metrics.ShannonMI_Zsh);
end

if iscell(R)
    for rr = 1:size(R,1)
        for cc = 1:size(R,2)
            if isempty(R{rr,cc})
                zMI(rr,cc) = NaN;
            else
                try
                    zMI(rr,cc) = nanmedian(R{rr,cc}.info_metrics.ShannonMI_Zsh);
                catch
                    try
                        zMI(rr,cc) = nanmedian(R{rr,cc}.info_metrics_1.ShannonMI_Zsh);
                    catch
                        zMI(rr,cc) = nanmedian(R{rr,cc}.MIs);
                    end
                end
            end
        end
    end
end
