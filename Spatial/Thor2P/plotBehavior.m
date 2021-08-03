function plotBehavior (ei)

m = matfile(ei.behaviorFile);
b = m.behavior;

nf = fields(b);

for ii = 1:length(nf)
    cmdTxt = sprintf('%s = b.%s;',nf{ii},nf{ii});
%     display(cmdTxt);
    eval(cmdTxt);
end

figure(1);clf;
subplot 211
plot(t,EA,t,EB); hold on;
plot(t,directionOfMotion + 2);
plot(t,reward + 4,'m');
ylim([-3 7])
subplot 212
% plot(tA,[0;speedA]);%,tA,[0;speedB]);
plot(t,thetaA);