function compare_spontaneous_activity

protocol_C = '10_C';
protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% ET_C = evalin('base',sprintf('ET%s',protocol_C));
% ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)
n = 0;
%%
out_C = get_spike_rate(ei_C);
out_A = get_spike_rate(ei_A);


n=0;
%%
if 1
    data_C = out_C.m_sp_animal_th;
    data_A = out_A.m_sp_animal_th;
    data = [data_C' data_A'];
    minBin = min([out_C.allVals_th;out_A.allVals_th]);
    maxBin = max([out_C.allVals_th;out_A.allVals_th])+20;
    incr = (maxBin-minBin)/50;
    hf = figure(1002);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = {'Ctrl','AD'};
    ylim([0 15]);
    xlim([0 maxBin+10])
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{length(legs)+1} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/3 dy/15];
    putLegend(ha,legs,'colors',colors);
    axes(ha);
    h = xlabel('Spike Rate');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0.1 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    file_name = fullfile(mData.pdf_folder,sprintf('%s_distribution spike rate',mfilename));
    save_pdf(hf,mData.pdf_folder,file_name,600);
end


n = 0;


function out = get_spike_rate(ei_C)
allVals = []; allVals_motion = []; allVals_rest = []; allVals_th = [];
for ii = 1:length(ei_C)
    ei = ei_C{ii};
    spSigAll = nan(length(ei.plane{1}.tP.deconv.caSigAll),length(ei.plane{1}.tP.deconv.caSigAll{1}));
    for cn = 1:length(ei.plane{1}.tP.deconv.caSigAll)
        spSigAll(cn,:) = ei.plane{1}.tP.deconv.spSigAll{cn}';
    end
    iscell = logical(ei.plane{1}.tP.iscell(:,1));
    spSigAll = spSigAll(iscell,:);
    speed = ei.b.fSpeed;
    inds_motion = find(speed > 0);
    inds_rest = find(speed == 0);
%     signal = spSigAll(1,:);
%     figure(100);clf;plot(signal);hold on;plot(ones(size(signal))*(mean(signal)+3*std(signal)))
    frames_f = ei.plane{1}.b.frames_f;
    [~,ind_frames_motion,~] = intersect(frames_f,inds_motion);
    [~,ind_frames_rest,~] = intersect(frames_f,inds_rest);
    m_sp_animal{ii} = nanmean(spSigAll,2);
    m_sp_animal_motion{ii} = nanmean(spSigAll(:,ind_frames_motion),2);
    m_sp_animal_rest{ii} = nanmean(spSigAll(:,ind_frames_rest),2);
    allVals = [allVals;m_sp_animal{ii}];
    allVals_motion = [allVals_motion;m_sp_animal_motion{ii}];
    allVals_rest = [allVals_rest;m_sp_animal_rest{ii}];
    m_sp_animal_level(ii) = nanmean(m_sp_animal{ii});
    m_sp_animal_level_motion(ii) = nanmean(m_sp_animal_motion{ii});
    m_sp_animal_level_rest(ii) = nanmean(m_sp_animal_rest{ii});
    thr = nanmean(spSigAll,2) + 3*nanstd(spSigAll,[],2);
    for cn = 1:size(spSigAll,1)
        m_sp_animal_th{ii}(cn,1) = nanmean(spSigAll(cn,spSigAll(cn,:) > thr(cn)));
    end
    m_sp_animal_level_th(ii) = nanmean(m_sp_animal_th{ii});
    allVals_th = [allVals;m_sp_animal_th{ii}];
end
out.m_sp_animal = m_sp_animal;
out.allVals = allVals;
out.m_sp_animal_level = m_sp_animal_level;
out.m_sp_animal_motion = m_sp_animal_motion;
out.allVals_motion = allVals_motion;
out.m_sp_animal_level_motion = m_sp_animal_level_motion;
out.m_sp_animal_rest = m_sp_animal_rest;
out.allVals_rest = allVals_rest;
out.m_sp_animal_level_rest = m_sp_animal_level_rest;

out.m_sp_animal_th = m_sp_animal_th;
out.m_sp_animal_level_th = m_sp_animal_level_th;
out.allVals_th = allVals_th;