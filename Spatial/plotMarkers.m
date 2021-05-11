function plotMarkers(b,markers1,markers2,fn,motion)

if ~exist('motion','var')
    motion = 0;
end

if motion
    plotMotion(b,markers1,markers2,fn);
    return;
end

figure(fn);clf;
% b = ei.b;

if isfield(b,'air_puff_raw')
    plot(b.ts,b.air_puff_raw,'color','b','linewidth',1.5);hold on;
end
hold on;
plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r','linewidth',1.5);
if isfield(b,'stim_raw')
    plot(b.ts,0.75*b.stim_raw+0,'color','g','linewidth',1.5);
end
plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m','linewidth',1.5);
xlabel('Time (sec)');
set(gca,'FontSize',12,'FontWeight','Bold');
if ~isempty(markers1)
    plot(b.ts(markers1),0.6*ones(size(markers1)),'ro');
    putText(b.ts(markers1),0.6*ones(size(markers1)))
end
if ~isempty(markers2)
    plot(b.ts(markers2),0.5*ones(size(markers2)),'mo');
    putText(b.ts(markers2),0.25*ones(size(markers2)))
end
if isempty(markers2)
    rightShift = 1;
    for ii = 1:length(markers1)
        text(b.ts(markers1(ii))+rightShift,0.9,num2str(ii));
    end
end
% xlim([b.ts(markers1(1))-15 b.ts(markers2(end))+15]);

function putText(xs,ys)
dx = xs(2)-xs(1);
for ii = 1:length(xs)
    text(xs(ii)+dx/3,ys(ii),num2str(ii));
end


function plotMotion(b,motionOnsets,motionOffsets,fn)
figure(fn);clf
plot(b.ts,b.fSpeed);hold on;
tMOn = zeros(size(b.fSpeed));
tMOn(motionOnsets) = 15;
plot(b.ts,tMOn,'r');
tMOn = zeros(size(b.fSpeed));
tMOn(motionOffsets) = 15;
plot(b.ts,tMOn,'k');
plot(b.ts,b.air_puff_raw,'k');
