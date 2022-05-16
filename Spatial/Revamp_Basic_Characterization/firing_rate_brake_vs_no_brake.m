function firing_rate_motion_vs_rest
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei_C = evalin('base','ei'); 
    for ii = 1:length(ei_C)
        [sig_B{ii},sig_NB{ii}] = get_spike_rate_max(ei_C{ii});
    end
    break
end
n = 0;
%%
sigB = exec_fun_on_cell_mat(sig_B,'mean');
sigNB = exec_fun_on_cell_mat(sig_NB,'mean');

[within,dvn,xlabels] = make_within_table({'Br'},[2]);
dataT = make_between_table({[sigB',sigNB']},dvn);
ra = RMA(dataT,within);
ra.ranova

%%
function [sig_B,sig_NB] = get_spike_rate(ei)
b = ei.b;

brakeOn1 = b.stim_r(1);
brakeOff1 = b.air_puff_f(10);

brakeOn2 = b.stim_r(21);
brakeOff2 = b.air_puff_f(50);

nbrakeOn = b.air_puff_r(11);
nbrakeOff = b.air_puff_f(40);

sig_B = [];
sig_NB = [];
for pp = 1:length(ei.plane)
    eip = ei.plane{pp};
    bp = eip.b;
    tspSigAll = [];
    this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
    tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);

    indsB1 = bp.frames_f > brakeOn1 & bp.frames_f < brakeOff1;
    indsB2 = bp.frames_f > brakeOn2 & bp.frames_f < brakeOff2;
    cellSigB = mean([mean(tspSigAll(:,indsB1),2) mean(tspSigAll(:,indsB2),2)],2);
    sig_B = [sig_B;cellSigB];
    
    indsB1 = bp.frames_f > nbrakeOn & bp.frames_f < nbrakeOff;
    cellSig_nB = mean(tspSigAll(:,indsB1),2);
    sig_NB = [sig_NB;cellSig_nB];
end
  
%%
function [sig_B,sig_NB] = get_spike_rate_max(ei)
b = ei.b;

brakeOn1 = b.stim_r(1);
brakeOff1 = b.air_puff_f(10);

brakeOn2 = b.stim_r(21);
brakeOff2 = b.air_puff_f(50);

nbrakeOn = b.air_puff_r(11);
nbrakeOff = b.air_puff_f(40);

sig_B = [];
sig_NB = [];
for pp = 1:length(ei.plane)
    eip = ei.plane{pp};
    bp = eip.b;
    tspSigAll = [];
    this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
    tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);

    indsB1 = bp.frames_f > brakeOn1 & bp.frames_f < brakeOff1;
    indsB2 = bp.frames_f > brakeOn2 & bp.frames_f < brakeOff2;
    cellSigB = max([max(tspSigAll(:,indsB1),[],2) max(tspSigAll(:,indsB2),[],2)],[],2);
    sig_B = [sig_B;cellSigB];
    
    indsB1 = bp.frames_f > nbrakeOn & bp.frames_f < nbrakeOff;
    cellSig_nB = max(tspSigAll(:,indsB1),[],2);
    sig_NB = [sig_NB;cellSig_nB];
end
  