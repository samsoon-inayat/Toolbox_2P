function processABF (ei)

if exist(ei.behaviorFile)
    display('Behavior file already exists');
    return;
end


abfFiles = dir(sprintf('%s\\*.abf',ei.dataFolder));
[abfd,abfsi,abfh]=abf2load(makeName(abfFiles(1).name,ei.dataFolder));
if size(abfd,2) == 5
    FT = abfd(:,1); % frame trigger
    EA = abfd(:,2); % encoder channel A
    EB = abfd(:,3); % encoder channel B
    ER = abfd(:,4); % encoder revolution
    reward = abfd(:,5); % reward
end
if size(abfd,2) == 6
    FT = abfd(:,1); % frame trigger
    ephys = abfd(:,2);
    EA = abfd(:,3); % encoder channel A
    EB = abfd(:,4); % encoder channel B
    ER = abfd(:,5); % encoder revolution
    reward = abfd(:,6); % reward
end

tA = linspace(0,abfsi*length(reward)*1e-6,length(reward))';

EA(EA<2.5) = 0;
EB(EB<2.5) = 0;
EA(EA>2.5) = 1;
EB(EB>2.5) = 1;
reward(reward < 5) = 0;
reward(reward > 5) = 1;
FT(FT<4.8) = 0;
FT(FT>0) = 1;

dEA = find(diff(EA) == -1);

indices  = EA == 0 & EB == 0;
bothZerosOffset = diff(indices) == -1;
bothZerosOnset = diff(indices) == 1;
ibzoff = find(bothZerosOffset == 1); 

directionOfMotion = zeros(size(EA)); % 0 for forward 1 for reverse if EB leads EA then its forward movement if EB lags EA then it is backward movement

if EA(ibzoff(1)+1)
    directionOfMotion(1:ibzoff(1)) = 1;
end

for ii = 1:(length(ibzoff)-1)
    if EA(ibzoff(ii)+1)
        directionOfMotion(ibzoff(ii):ibzoff(ii+1)) = 1;
    end
    if EB(ibzoff(ii)+1)
        directionOfMotion(ibzoff(ii):ibzoff(ii+1)) = 0;
    end
end

rewardIndices = find(diff(reward)==-1);
framesIndices = find(diff(FT) == 1);

EBOnset = find(diff(EB) == 1);
angularPosition_pulseCountB = zeros(size(EA));
currentCount = 1;
for ii = 1:(length(EBOnset)-1)
    if directionOfMotion(EBOnset(ii)) == 1
        currentCount = currentCount - 1;
        angularPosition_pulseCountB(EBOnset(ii):EBOnset(ii+1)) = currentCount;
        continue;
    end
    if directionOfMotion(EBOnset(ii)) == 0
        angularPosition_pulseCountB(EBOnset(ii):EBOnset(ii+1)) = currentCount;
        currentCount = currentCount+1;
    end
end
angularPosition_pulseCountB(EBOnset(end):length(EB)) = currentCount;
thetaB = angularPosition_pulseCountB*2*pi/500;

EAOnset = find(diff(EA) == 1);
angularPosition_pulseCountA = zeros(size(EA));
currentCount = 1;
for ii = 1:(length(EAOnset)-1)
    if directionOfMotion(EAOnset(ii)) == 1
        currentCount = currentCount - 1;
        angularPosition_pulseCountA(EAOnset(ii):EAOnset(ii+1)) = currentCount;
        continue;
    end
    if directionOfMotion(EAOnset(ii)) == 0
        angularPosition_pulseCountA(EAOnset(ii):EAOnset(ii+1)) = currentCount;
        currentCount = currentCount+1;
    end
end
angularPosition_pulseCountA(EAOnset(end):length(EA)) = currentCount;
thetaA = angularPosition_pulseCountA*2*pi/500;

% linearDistanceA = 5 * thetaA;
% linearDistanceB = 5 * thetaB;

% dT = diff(tA);
% dThetaA = diff(thetaA);
% wA = dThetaA./dT;
% % speedA  = 5*wA;
% dThetaB = diff(thetaB);
% wB = dThetaB./dT;
% speedB = 5*wB;

behavior.t = tA;
behavior.directionOfMotion = directionOfMotion;
behavior.thetaA = thetaA;
behavior.thetaB = thetaB;
behavior.rewards = tA(rewardIndices);
behavior.rewardIndices = rewardIndices;
behavior.frames = tA(framesIndices);
behavior.framesIndices = framesIndices;
behavior.EA = EA;
behavior.EB = EB;
behavior.reward = reward;
if size(abfd,2) == 6
    behavior.ephys = ephys;
end

m = matfile(ei.behaviorFile);
m.behavior = behavior;


% figure(1);clf;
% subplot 211
% plot(tA,EA,tA,EB); hold on;
% plot(tA,directionOfMotion + 2);
% plot(tA,reward + 4,'m');
% ylim([-3 7])
% subplot 212
% % plot(tA,[0;speedA]);%,tA,[0;speedB]);
% plot(tA,thetaA);

