function save_data_for_NB_Decoding

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei'); 


selContexts = [3 4 5];
rasterNames = {'airT','airT','airT'};
raster_data_C = get_rasters_data(ei_C,selContexts,rasterNames);

selContexts = [3 4 5];
rasterNames = {'airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
% view_population_vector(Rs,mRs,300);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(RsC);

pcs = 2;

for an = 1:5
    for cn = 1:3
        if pcs == 2 % not place cells
            [DistC{an,cn},TC{an,cn},SpaceC{an,cn},RastersC{an,cn},SpeedC{an,cn},frame_rateC(an,1)] = getXYs1(raster_data_C{an,cn},~resp_valsC{an}(:,cn));
        end
        if pcs == 1 % place cells
            [DistC{an,cn},TC{an,cn},SpaceC{an,cn},RastersC{an,cn},SpeedC{an,cn},frame_rateC(an,1)] = getXYs1(raster_data_C{an,cn},resp_valsC{an}(:,cn));
        end
        if pcs == 0 % all cells
            [DistC{an,cn},TC{an,cn},SpaceC{an,cn},RastersC{an,cn},SpeedC{an,cn},frame_rateC(an,1)] = getXYs1(raster_data_C{an,cn},[]);
        end
    end
end

if pcs == 2
    fileName = fullfile(mData.pd_folder,sprintf('NB_decoding_data_not_place_cells.mat'));
end
if pcs == 1
    fileName = fullfile(mData.pd_folder,sprintf('NB_decoding_data_place_cells.mat'));
end
if pcs == 0
    fileName = fullfile(mData.pd_folder,sprintf('NB_decoding_data_all_cells.mat'));
end
save(fileName,'DistC','TC','SpaceC','RastersC','SpeedC','frame_rateC','-v7.3');
return



train = 1:10;% test = 2:2:10;
% 
% [outNB.aXs_C_train,outNB.aYs_C_train] = getXYs1(raster_data_C,cn,train,pMs_C); [outNB.aXs_C_test,outNB.aYs_C_test] = getXYs1(raster_data_C,cn,test,pMs_C);
% 
% [outNB.aXs_A_train,outNB.aYs_A_train] = getXYs1(raster_data_A,cn,train,pMs_A); [outNB.aXs_A_test,outNB.aYs_A_test] = getXYs1(raster_data_A,cn,test,pMs_A);

for cn = 1:4
    [outNB.XsC{cn},outNB.YsC{cn}] = getXYs1(raster_data_C,cn,train,pMs_C);
    [outNB.XsA{cn},outNB.YsA{cn}] = getXYs1(raster_data_A,cn,train,pMs_A);
end

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

function [dist,t,space,response,speed,frame_rate] = getXYs1(raster_data,resp)
rd = raster_data.fromFrames;
dist = rd.dist;
t = rd.duration;
space = rd.space;
if isempty(resp)
    response = rd.sp_rasters;
else
    response = rd.sp_rasters(:,:,logical(resp));
end
frame_rate = raster_data.thorexp.frameRate;
speed = rd.speed;
n = 0;





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
