function findResponsivenessOfTrials (ags,pags,stim_folders_info,indices,frames_info,maskFactor)

for ii = 1:length(indices)
    indx = indices(ii)
    preProcessFolder(ags,pags,stim_folders_info,indx,frames_info,maskFactor);
end


function preProcessFolder(ags,pags,stim_folders_info,indx,frames_info,maskFactor)
gg = stim_folders_info.list(indx,2);
aa = stim_folders_info.list(indx,3);
rr = stim_folders_info.list(indx,4);
thisRecording = ags(gg).animals{aa}.eRecordings{rr};
pThisRecording = pags(gg).animals{aa}.eRecordings{rr};
disp(thisRecording.root_folder);
disp(pThisRecording.root_folder);
mask = getMask(ags(gg).animals{aa}.root_folder,maskFactor);
[rois,roi_names] = getROIs_VSDI(ags,gg,aa);
% fileName = makeName('rois.mat',pags(gg).animals{aa}.root_folder);
% save(fileName,'rois','-v7.3');

hitrials = thisRecording.hiTrials;
lotrials = thisRecording.loTrials;
notrials = thisRecording.noTrials;
disp('Loading data');

trials = notrials; noallImg = cell(1,length(trials));
root_folder = pThisRecording.root_folder;
for ii = 1:length(trials)
    filenamedfSeq = makeName([trials{ii} '.mat'],root_folder);
    temp = load(filenamedfSeq);
    noallImg{ii}.img = temp.df;
end
trials = hitrials; hiallImg = cell(1,length(trials));
for ii = 1:length(trials)
    filenamedfSeq = makeName([trials{ii} '.mat'],root_folder);
    temp = load(filenamedfSeq);
    hiallImg{ii}.img = temp.df;
end
trials = lotrials; loallImg = cell(1,length(trials));
for ii = 1:length(trials)
    filenamedfSeq = makeName([trials{ii} '.mat'],root_folder);
    temp = load(filenamedfSeq);
    loallImg{ii}.img = temp.df;
end

disp('Processing data');
sel_rois = 1:19;
for tt = 1:length(trials)
    for rrr = 1:length(sel_rois)
        rr = sel_rois(rrr);
%         FRVhi(rr,tt) = findFRV1(hiallImg{tt}.img(rr,:),noallImg{tt}.img(rr,:),frames_info);
%         FRVlo(rr,tt) = findFRV1(loallImg{tt}.img(rr,:),noallImg{tt}.img(rr,:),frames_info);
%         pause;
%         [FRVhi(rrr,tt),~,~] = findFRV(hiallImg{tt}.img(rr,:),frames_info);
%         [FRVlo(rrr,tt),~,~] = findFRV(loallImg{tt}.img(rr,:),frames_info);
%         [FRVno(rrr,tt),~,~] = findFRV(noallImg{tt}.img(rr,:),frames_info);
        resp.FRVhijb(rrr,tt) = jbtest(hiallImg{tt}.img(rr,:));
        resp.FRVlojb(rrr,tt) = jbtest(loallImg{tt}.img(rr,:));
        resp.FRVnojb(rrr,tt) = jbtest(noallImg{tt}.img(rr,:));
        resp.FRVhili(rrr,tt) = lillietest(hiallImg{tt}.img(rr,:));
        resp.FRVloli(rrr,tt) = lillietest(loallImg{tt}.img(rr,:));
        resp.FRVnoli(rrr,tt) = lillietest(noallImg{tt}.img(rr,:));
%         resp.FRVhiad(rrr,tt) = adtest(hiallImg{tt}.img(rr,:));
%         resp.FRVload(rrr,tt) = adtest(loallImg{tt}.img(rr,:));
%         resp.FRVnoad(rrr,tt) = adtest(noallImg{tt}.img(rr,:));
    end
end
fileName = makeName('trials_responsiveness.mat',root_folder);
save(fileName,'-struct','resp');

% mFRVno = mean(FRVno,2);
% stdFRVno = std(FRVno,[],2);
% quit = 0
% fact = 0;
% resphi = zeros(length(sel_rois),length(trials));
% resplo = zeros(length(sel_rois),length(trials));
% for tt = 1:length(trials)
%     for rr = 1:length(sel_rois)
%         if FRVhi(rr,tt) > FRVno(rr,tt)% + fact * stdFRVno(rr)
%             resphi(rr,tt) = 1;
%         end
%         if FRVlo(rr,tt) > FRVno(rr,tt)% + fact * stdFRVno(rr)
%             resplo(rr,tt) = 1;
%         end
%     end
% end
% 
% for rrr = 1:length(sel_rois)
%     rr = sel_rois(rrr);
%     for tt = 1:length(trials)
%         if quit
%             break;
%         end
%         figure(1);clf;
%         subplot 211
%         plot(hiallImg{tt}.img(rr,:));hold on;
%         plot(noallImg{tt}.img(rr,:),'r');hold on;
%         titleText = sprintf('%s -- Trial - %d -- Resp - %d -- FRV [%.2f %.2f]',roi_names{rr},tt,resphi(rrr,tt),FRVhi(rrr,tt),FRVno(rrr,tt));
%         title(titleText);
%         subplot 212
%         plot(loallImg{tt}.img(rr,:));hold on;
%         plot(noallImg{tt}.img(rr,:),'r');hold on;
%         titleText = sprintf('%s -- Trial - %d -- Resp - %d -- FRV [%.2f %.2f]',roi_names{rr},tt,resplo(rrr,tt),FRVlo(rrr,tt),FRVno(rrr,tt));
%         title(titleText);
%         ch = getkey;
%         if ch == 27
%             quit = 1;
%             break;
%         end
%     end
% end

disp('Done with this folder');

function FRV = findFRV1(sig,signo,frames_info)
hsig = jbtest(sig);
hsigno = jbtest(signo);
if hsig & ~hsigno
    FRV = 1;
else
    FRV = 0;
end
return;
stimFrame = find(frames_info(4):frames_info(5) == frames_info(1))+2;
st1 = sig(1:stimFrame);
st2 = sig((stimFrame+1):end);
means1 = mean(st1)*ones(size(st1));
err1 = std(st1)*ones(size(st1));
err2 = std(st2)*ones(size(st2));
means2 = mean(st2)*ones(size(st2));
means = [means1 means2];
err = [err1 err2];
xs = 1:length(means);

st1 = signo(1:stimFrame);
st2 = signo((stimFrame+1):end);
means1 = mean(st1)*ones(size(st1));
means2 = mean(st2)*ones(size(st2));
meansno = [means1 means2];
err1 = std(st1)*ones(size(st1));
err2 = std(st2)*ones(size(st2));
errno = [err1 err2];
xsno = 1:length(meansno);

figure(1);clf;
plot(sig);hold on;
y = means;
x = xs;
errr = err;
plot(x,y,'c');
patch([x fliplr(x)],[y+errr fliplr(y-errr)],'c','FaceAlpha',0.05,'EdgeColor','None');
plot(signo,'r');
y = meansno;
x = xsno;
errr = errno;
plot(x,y,'m');
% patch([x fliplr(x)],[y+errr fliplr(y-errr)],'m','FaceAlpha',0.05,'EdgeColor','None');
[FRV ft1 ft2] = findFRV(sig,frames_info);
[FRVno ft1no ft2no] = findFRV(signo,frames_info);
titleText = sprintf('FRV [%.2f %.2f %.2f] - FRV no [%.2f %.2f %.2f]',FRV,ft1,ft2,FRVno,ft1no,ft2no);
title(titleText);

FRV = 1;



function [FRV ft1 ft2] = findFRV(sig,frames_info)
extra = 2;
stimFrame = find(frames_info(4):frames_info(5) == frames_info(1))+extra;
theEnd = length(sig);%(stimFrame+1)+15;
st1 = sig(1:stimFrame);
st2 = sig((stimFrame+1):theEnd);
ft1 = var(st1)/abs(mean(st1));
ft2 = var(st2)/abs(mean(st2));
fall = var(sig(1:theEnd))/abs(mean(sig(1:theEnd)));
FRV = (ft2)/ft1;
disp([ft1 ft2 FRV])




