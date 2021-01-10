function save_data_for_NB_Decoding

protocol_C = '10_C';
protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C);
selAnimals_A = 1:length(ei_A);

[out.aXs_C,out.aYs_C] = getXYs(ei_C,selAnimals_C);
[out.aXs_A,out.aYs_A] = getXYs(ei_A,selAnimals_A);
save('NB_decoding.mat','-struct','out','-v7.3');
n = 0;

function [aXs,aYs] = getXYs(ei_C,selAnimals_C)

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
    end
    aXs{ii} = Xs;
    aYs{ii} = Ys;
end
