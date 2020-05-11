function st_demo_ImagingComp
ei = evalin('base','ei{1}');

% Demo program
% STeP is applied to simulated data
%
% Copyright (C) 2016, Yusuke Takeda, ATR, takeda@atr.jp

% Set parameters for this simulation test
% 
% clear all
% close all
load('coi.mat');
signals = ei.tP.signals;

data = signals(ei.areCells(coi),:)';


data_length=20000;% Length of simulated data
N=50;% Length of spatiotemporal pattern
Npattern=10;% Number of spatiotemporal patterns

%% Estimate spatiotemporal patterns and their onsets using STeP

% Apply STeP to simulated data
% [eonset,epattern]=st_STeP(data,N,Npattern);
load('STeP_eonset_epattern.mat');
% Show estimation result
figure
eonset_timeseries=st_make_onset_timeseries(eonset,data_length);
for pt=1:Npattern
    subplot(Npattern,5,5*(pt-1)+1)
    imagesc(squeeze(epattern(:,pt,:))')
    colorbar
    if pt==1
        title('Estimated spatiotemporal pattern')
    end
    if pt==Npattern
        xlabel('Time')
    end
    ylabel('Channel')
    subplot(Npattern,5,5*(pt-1)+2:5*pt)
    plot(eonset_timeseries(:,pt))
    if pt==1
        title('Estimated onset timeseries')
    end
    if pt==Npattern
        xlabel('Time')
    end
end

% Adjust estimation result 
% (reorder estimated spatiotemporal patterns and shift estimated onsets along time)
[aonset,apattern]=st_adjust(data,pattern,eonset,N);

% Show adjusted result
figure
aonset_timeseries=st_make_onset_timeseries(aonset,data_length);
for pt=1:Npattern
    subplot(Npattern,5,5*(pt-1)+1)
    imagesc(squeeze(apattern(:,pt,:))')
    colorbar
    if pt==1
        title('Adjusted spatiotemporal pattern')
    end
    if pt==Npattern
        xlabel('Time')
    end
    ylabel('Channel')
    subplot(Npattern,5,5*(pt-1)+2:5*pt)
    plot(aonset_timeseries(:,pt))
    if pt==1
        title('Adjusted onset timeseries')
    end
    if pt==Npattern
        xlabel('Time')
    end
end

%% Quantify estimation accuracy

% Correlation coefficient between true and estimated spatiotemporal patterns
r1=zeros(Npattern,1);
for pt=1:Npattern
    true=pattern(:,pt,:);
    estimated=apattern(:,pt,:);
    r1(pt,1)=st_cc(true(:),estimated(:));
end
r=mean(r1);
fprintf('Correlation coefficient is %1.2f.\n',r)

% Normalized distance from true onsets
d1=zeros(Nonset,1);
nd1=zeros(Npattern,1);
for pt=1:Npattern
    aonset1=aonset(aonset(:,pt)>0,pt);
    mean_ioi=mean(diff(sort(aonset1)));
    for on=1:Nonset
        d1(on,1)=min(abs(onset(on,pt)-aonset(:,pt)));
    end
    nd1(pt,1)=mean(d1)/mean_ioi;
end
nd=mean(nd1);
fprintf('Normalized distance from true onset is %1.4f.\n',nd)

% Normalized number of estimated onsets
nn1=zeros(Npattern,1);
for pt=1:Npattern
    a=find(aonset(:,pt)>0);
    nn1(pt,1)=length(a)/Nonset;
end
nn=mean(nn1);
fprintf('Normalized number of estimated onsets is %1.2f.\n',nn)

