function pol = find_resp_polarity(Rs,good_FR)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        resp = good_FR{rr,cc};
        if strcmp(R.marker_name,'air55T')
            tpol = find_pol(R,resp);
        else
            disp('some problem');
            error;
        end
        
    end
end

function tpol = find_pol(R,resp)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
markerType = R.marker_name;
timeBefore = str2num(markerType(end-1));
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
column_index = round(timeBefore * SR);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;

% resp = (R.info_metrics.ShannonMI_Zsh > 1.96)';
% return;

group = [];
for ii = 1:size(cis,2)
    group = [group ii*ones(1,length(cis(1,ii):cis(2,ii)))];
end
p = NaN(size(rasters,3),1);
hv = p;
resp = logical(zeros(size(p)));
excinh = p;
for ii = 1:size(rasters,3)
    thisRaster = rasters(:,:,ii);
    m_thisRaster = nanmean(thisRaster);
    [p(ii),atabk,statsk] = kruskalwallis(m_thisRaster,group,'nodisplay');
    vert = nansum(thisRaster,2);
    hv(ii) = sum(vert>0) > 5;
    if p(ii) < 0.05 & hv(ii) == 1
        resp(ii) = 1;
        if mean(m_thisRaster(group==1)) > mean(m_thisRaster(group==2)) && mean(m_thisRaster(group==3)) > mean(m_thisRaster(group==2))
            excinh(ii) = 0;
        else
            excinh(ii) = 1;
        end
    end
end

resp = resp';

