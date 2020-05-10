function behaviorPlots(ei)

ei = evalin('base','ei{1}');

b = ei.b;
names = fieldnames(b);
for ii = 1:length(names)
    eval(sprintf('%s = b.%s;', names{ii},names{ii}));
end
clear b;

% % photo_sensor_f = validate_photo_sensor_signal(photo_sensor_f);
for ii = 1:length(trials)
    st = air_puff_r(ii);
    se = air_puff_f(ii);
    trialSpeeds = speed(st:se);
    figure(300);clf;
    plot(trialSpeeds);
    n = 0;
    TS(ii) = mean(trialSpeeds);
    eTS(ii) = std(trialSpeeds)/sqrt(length(trialSpeeds));
    trialDuration(ii) = ts(se)-ts(st);
    aTS(ii) = (dist(se)-dist(st))/(ts(se)-ts(st));
end

for ii = 1:(length(trials)-1)
    st = air_puff_f(ii);
    se = air_puff_r(ii+1);
    trialSpeeds = speed(st:se);
    figure(300);clf;
    plot(trialSpeeds);
    nTS(ii) = mean(trialSpeeds);
    enTS(ii) = std(trialSpeeds)/sqrt(length(trialSpeeds));
end

figure(300);clf;
subplot 211
plot(aTS,'o-','linewidth',2);hold on;
% errorbar(TS,eTS);
xlim([0 length(trials)+1]);
xlabel('Trials');
ylabel('Average Speed (cm/sec)')
set(gca,'Linewidth',2,'box','off','TickDir','out','FontSize',11,'FontWeight','Bold')

% plot(nTS);
% errorbar(nTS,enTS);
subplot 212
plot(trialDuration,'o-','linewidth',2);
xlabel('Trials');
ylabel('Trial Duration (sec)');
set(gca,'Linewidth',2,'box','off','TickDir','out','FontSize',11,'FontWeight','Bold')
xlim([0 length(trials)+1]);
% figure(300);clf
% plot(ts,speed,'r');hold on;
% % plot(ts,air_puff_raw*max(speed),'b','linewidth',2);
% dsDist = diff([dist]);
% dts = diff([ts]);
% eSpeed = dsDist./dts;
% plot(ts(2:end)+10,max(speed)*dsDist/max(dsDist),'k');
n=0;
if ~isfield(ei.db,'brain_location')
    pdfFileName = sprintf('behavior_%s_%s_%s.pdf',ei.animal_id,'RSC',ei.exp_date);
else
    pdfFileName = sprintf('behavior_%s_%s_%s.pdf',ei.animal_id,ei.db.brain_location,ei.exp_date);
end
save2pdf(pdfFileName,gcf,600);

function photoS = validate_photo_sensor_signal(photoS)
idx = find(diff(photoS)<10000);
photoS(idx+1) = [];
