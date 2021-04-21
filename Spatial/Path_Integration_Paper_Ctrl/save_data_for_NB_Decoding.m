function save_data_for_NB_Decoding

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;

% typeP = 'Population';cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN; FR = NaN;
typeP = 'Place';cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [1 150]; rs_th = 0.3; FR = NaN;%[0.1 5000];

conditionsAndRasterTypes = [11;21;31;41];
% conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out = read_data_from_base_workspace_onlyC(selC)

ei_C = out.eis{1}; %ei_A = out.eis{2};
pMs_C = out.pMs{1}; %pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; %paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; %cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; %selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; %perc_cells_A = out.perc_cells{2};

% [out.aXs_C,out.aYs_C,out.aYs1_C] = getXYs(ei_C,selAnimals_C);
% [out.aXs_A,out.aYs_A,out.aYs1_A] = getXYs(ei_A,selAnimals_A);
% 
% fileName = fullfile(mData.pd_folder,sprintf('NB_decoding.mat'));
% save(fileName,'-struct','out','-v7.3');
% return;
selContexts = [1 2 3 4];
% rasterNames = {'airD','airD','airD','airD'};
rasterNames = {'airT','airT','airT','airT'};
raster_data_C = get_rasters_data(ei_C,selContexts,rasterNames);
% raster_data_A = get_rasters_data(ei_A,selContexts,rasterNames);

cn = 1; train = 1:2:10; test = 2:2:10;
% 
[outNB.aXs_C_train,outNB.aYs_C_train] = getXYs1(raster_data_C,cn,train,pMs_C); [outNB.aXs_C_test,outNB.aYs_C_test] = getXYs1(raster_data_C,cn,test,pMs_C);

% [outNB.aXs_A_train,outNB.aYs_A_train] = getXYs1(raster_data_A,cn,train,pMs_A); [outNB.aXs_A_test,outNB.aYs_A_test] = getXYs1(raster_data_A,cn,test,pMs_A);

fileName = fullfile(mData.pd_folder,sprintf('NB_decoding.mat'));
save(fileName,'-struct','outNB','-v7.3');
return;

fileName = fullfile(mData.pd_folder,sprintf('NB_decoding_C_%s.mat',typeP));
save(fileName,'-struct','all_out_C_or','-v7.3');
fileName = fullfile(mData.pd_folder,sprintf('NB_decoding_A_%s.mat',typeP));
save(fileName,'-struct','all_out_A_or','-v7.3');
n = 0;

n = 0;
% filename_or = fullfile(mData.pd_folder,sprintf('pop_corr_mat_remapping_or_AD_%s.mat',typeP));
% if 0
%     cellSel_C = ''; cellSel_A = ''; cellSel_C = cpMs_C; cellSel_A = cpMs_A;
%     a_trials = {3:10};
%     for ii = 1:length(a_trials)
%         out = get_mean_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,cellSel_C,a_trials{ii});
%         [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an,out.all_corr_cell_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,min(out.sz(:)),cellSel_C);
%         all_out_C{ii} = out;
%         out = get_mean_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,cellSel_A,a_trials{ii});
%         [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an,out.all_corr_cell_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,49,cellSel_A);
%         all_out_A{ii} = out;
%     end
%     save(filename_or,'all_out_C','all_out_A','a_trials','-v7.3');
%     disp('Done');
%     return;
% else
%     temp = load(filename_or);
%     all_out_C_or = temp.all_out_C{1};     all_out_A_or = temp.all_out_A{1}; a_trials = temp.a_trials{1};
% end
% 

function [aXs_train,aYs_train] = getXYs1(raster_data,cn,train,pMs)
pMsC = pMs{cn};
for ii = 1:length(raster_data)
    iscell = raster_data{ii}{cn}.iscell;% & pMsC.cellSel{ii};
    rd = raster_data{ii}{cn}.fromFrames;
    
    Ys = [];
    dist = rd.dist(train,:);
    dist = reshape(dist',size(dist,1)*size(dist,2),1);
    inds = ~isnan(dist);
    Ys(:,1) = dist(inds);
   
    dist = rd.speed(train,:);
    dist = reshape(dist',size(dist,1)*size(dist,2),1);
    Ys(:,2) = dist(inds);
    if sum(isnan(Ys(:)))>0
        n = 0;
    end
    [rnanYs,cnanYs] = find(isnan(Ys));
    aYs_train{ii} = Ys';
    
    sp_rasters = permute(rd.sp_rasters(train,:,:),[2 1 3]);
    sp_rastersR = reshape(sp_rasters,size(sp_rasters,1)*size(sp_rasters,2),size(sp_rasters,3));
    sp_rastersR = sp_rastersR(:,iscell==1);
    
    
    Xs = sp_rastersR(inds,:);
    [rnanXs,cnanXs] = find(isnan(Xs));
    Xs(:,cnanXs) = [];
    Xs = fillmissing(Xs,'linear','EndValues','nearest');
%     mXs = nanmean(Xs,2);
    [Xs,mu,sigma] = zscore(Xs,0,2);
    if sum(isnan(Xs(:))) > 0
        n = 0;
    end
%     disp(max(Xs(:)));
%     [rs,cs] = find(Xs>1000);
%     if ~isempty(rs)
%          n = 0;
%     end
    aXs_train{ii} = Xs';
   
end

function [aXs,aYs,aYs1] = getXYs(ei_C,selAnimals_C)

for ii = 1:length(selAnimals_C)
    tei = ei_C{selAnimals_C(ii)};
    b = tei.b;
    a_dist = b.dist;
    a_dist1 = nan(size(a_dist));
    ps = b.photo_sensor_f;
    d_adist = diff(a_dist(ps));
    lengs =[];
    for psii = 2:length(ps)
        tps = ps(psii);
        lengs(psii-1) = a_dist(ps(psii))-a_dist(ps(psii-1));
    end
    for psii = 1:length(ps)
        if psii == 1
            if a_dist(ps(psii)) < max(lengs)
                a_dist1(1:(ps(psii)-1)) = a_dist(1:(ps(psii)-1)) + (round(median(lengs)) - a_dist((ps(psii)-1)));
                a_dist1(ps(psii)) = round(median(lengs));
                continue;
            end
        else
            a_dist1((ps(psii-1)+1):(ps(psii))) = a_dist((ps(psii-1)+1):(ps(psii)))-a_dist((ps(psii-1)));
        end
    end
    
    for pp = 1%:length(tei.plane)
        iscell = tei.plane{pp}.tP.iscell;
        frames_f = tei.plane{pp}.b.frames_f;
        frames_f = frames_f(find(frames_f < ps(end)));
        spSigAll = nan(length(tei.plane{1}.tP.deconv.spSigAll),length(frames_f));
        spSig = tei.plane{pp}.tP.deconv.spSigAll;
        for cc = 1:length(spSig)
            thisSig = spSig{cc}';
            spSigAll(cc,:) = thisSig(1:length(frames_f));
        end
        Xs{pp} = spSigAll(find(iscell(:,1)),:);
        Ys{pp} = a_dist1(frames_f);
        Ys1{pp} = b.air_puff_raw(frames_f);
    end
    aXs{ii} = Xs;
    aYs{ii} = Ys;
    aYs1{ii} = Ys1;
end
