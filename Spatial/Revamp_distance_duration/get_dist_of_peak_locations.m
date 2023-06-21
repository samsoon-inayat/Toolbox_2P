function [dist, distL, distTr, distTrL, allMI, allMIL, allMITr, allMITrL] = get_dist_of_peak_locations(cellpop,si,mRs,allpeakL_trials,allresp_trials,binEs,propsG)

for ii = 1:size(mRs,2)
    size_raster_2(1,ii) = size(mRs{1,ii},2);
end

if strcmp(cellpop,'resp')
    rccc = 1;
    trials = 1:10;
end

if strcmp(cellpop,'conj')
    rccc = 2;
    trials = 1:9;
end

if strcmp(cellpop,'comp1')
    rccc = 3;
    trials = 1:9;
end


if strcmp(cellpop,'comp2')
    rccc = 4;
    trials = 2:10;
end

n = 0;
%% plot average distributions of the peak locations themselves
for cn = 1:length(si)
    for an = 1:5
        pLs = allpeakL_trials{an,cn}; resp = allresp_trials{an,cn}; MIs = propsG.MI_trials{an,cn};
        trNi = 1;
        for trN = trials
            switch rccc
                case 1
                    t_resp = resp(:,trN);
                case 2
                    t_resp = resp(:,trN) & resp(:,trN+1);
                case 3
                    t_resp = resp(:,trN) & ~resp(:,trN+1);
                case 4
                    t_resp = ~resp(:,trN-1) & resp(:,trN);
            end
            pL_vals = pLs(:,trN); tMIs = MIs(t_resp,trN); 
            pL_vals = pL_vals(t_resp)/size_raster_2(cn);
            [bar1,binEs,binVals] = histcounts(pL_vals,binEs,'Normalization','probability');
            dist(cn,an,trNi,:) = bar1;
            if sum(bar1==0) > 0
                n = 0;
            end
            allMI(cn,an,trNi,:) = accumarray(binVals,tMIs,size(bar1'),[],0)'./accumarray(binVals,1,size(bar1'),[],0)';
%             err
            trNi = trNi + 1;
%             pL_vals_trials = [pL_vals_trials;bar1]; pL_vals_trialsLin = [pL_vals_trialsLin bar1];
        end
    end
end


if strcmp(cellpop,'resp')
    distT = permute(dist,[4 3 1 2]); % pL,tr,cn,an
    distL = (reshape(distT,500,5))';
    distTr = squeeze(mean(distT,2));
    distTrL = (reshape(distTr,50,5))';
    
    distT = permute(allMI,[4 3 1 2]); % pL,tr,cn,an
    allMIL = (reshape(distT,500,5))';
    allMITr = squeeze(mean(distT,2));
    allMITrL = (reshape(allMITr,50,5))';
end

if strcmp(cellpop,'conj')
    distT = permute(dist,[4 3 1 2]); % pL,tr,cn,an
    distL = (reshape(distT,450,5))';
    distTr = squeeze(mean(distT,2));
    distTrL = (reshape(distTr,50,5))';
    
    distT = permute(allMI,[4 3 1 2]); % pL,tr,cn,an
    allMIL = (reshape(distT,450,5))';
    allMITr = squeeze(mean(distT,2));
    allMITrL = (reshape(allMITr,50,5))';
end

if strcmp(cellpop,'comp1')
    distT = permute(dist,[4 3 1 2]); % pL,tr,cn,an
    distL = (reshape(distT,450,5))';
    distTr = squeeze(mean(distT,2));
    distTrL = (reshape(distTr,50,5))';
    
    distT = permute(allMI,[4 3 1 2]); % pL,tr,cn,an
    allMIL = (reshape(distT,450,5))';
    allMITr = squeeze(mean(distT,2));
    allMITrL = (reshape(allMITr,50,5))';
end

if strcmp(cellpop,'comp2')
    distT = permute(dist,[4 3 1 2]); % pL,tr,cn,an
    distL = (reshape(distT,450,5))';
    distTr = squeeze(mean(distT,2));
    distTrL = (reshape(distTr,50,5))';
    
    distT = permute(allMI,[4 3 1 2]); % pL,tr,cn,an
    allMIL = (reshape(distT,450,5))';
    allMITr = squeeze(mean(distT,2));
    allMITrL = (reshape(allMITr,50,5))';
end


