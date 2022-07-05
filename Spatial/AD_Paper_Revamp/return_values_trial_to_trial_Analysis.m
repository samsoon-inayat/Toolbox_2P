function [allRsC,allmRsT,allresp,pcs,event_type] = return_values_trial_to_trial_Analysis(o,si)
%%

%% find spatial trial to trial correlation
    trialNums = [1:10];
%    si = [C1_t_D C1_i_T C2_t_D C2_i_T C3_t_D C3_i_T C4_t_D C4_i_T];
   ind = 1;
   for rr = 1:4
       for cc = 1:10
           event_type{ind} = sprintf('D%d%d',rr,cc-1);
           ind = ind + 1;
       end
       for cc = 1:10
           event_type{ind} = sprintf('T%d%d',rr,cc-1);
           ind = ind + 1;
       end
   end
   
    Rs = o.Rs(:,si);mR = o.mR(:,si); RsG = Rs; siG = si; propsG = get_props_Rs(RsG,[40,100]); respG = propsG.vals;
    avgProps = get_props_Rs(Rs,[40,100]); respM = avgProps.good_FR;
    for cn = 1:length(si)
        trials = mat2cell([1:10]',ones(size([1:10]')));
        trials = mat2cell([trialNums]',ones(size([trialNums]')));
        RsC = repmat(Rs(:,cn),1,10);
        mRsCT = cell(size(RsC,1),length(trials));
        for ii = 1:length(trials)
            ii;
            [mRsCT(:,ii),~] = calc_mean_rasters(RsC(:,1),trials{ii});
        end
        allmRsT{cn} = mRsCT;
        allRsC{cn} = RsC;
    end
    disp('Done');
    %%



%% Overlap Indices ImageSC all

    avgProps = get_props_Rs(RsG,[50,100]); 
    respG = avgProps.vals;
    an  = 1:5; eic = 1; sp = 0; intersect_with_global = 0; only_global = 0;
    allresp = []; ind = 1;
    all_peakL = [];
    for cn = 1:length(si)
        mRsCT = allmRsT{cn};
        resp = []; peak_locations = [];
        for rr = 1:size(mRsCT,1)
            for cc = 1:size(mRsCT,2)
                this_mat = mRsCT{rr,cc};
                [~,peakL] = max(this_mat,[],2);
%                 size_tmat(rr,cc) = size(this_mat,2);
                resp{rr,cc} = sum(this_mat,2) > 0;
                if intersect_with_global
                    resp{rr,cc} = resp{rr,cc} & respG{rr,cn};
                end
                if only_global
                    resp{rr,cc} = respG{rr,cn};
                end
                if sp == 1
                    if cn == 1
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                    end
                    if cn == 2
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                    end
                    if cn == 3
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,3};
                    end
                    if cn == 4
                        respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                    end
                    if cn == 5
                        respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                    end
                end
                peakL(~resp{rr,cc}) = NaN;
                peak_locations{rr,cc} = peakL;
                if rr == 1
                    txl{ind} = sprintf('C%dT%d',cn,cc);
                    ind = ind + 1;
                end
            end
%             oc(rr,cn) = find_cells_based_on_cluster(cell2mat(resp(rr,:)));
        end
        allresp = [allresp resp]; all_peakL = [all_peakL peak_locations];
        
    end
    i_allresp = cell_list_op(allresp,[],'not');

    allrespOR = cell_list_op(allresp,[],'or',1);
    allrespAND = cell_list_op(allresp,[],'and',1);
    
    pallrespOR = 100*exec_fun_on_cell_mat(allrespOR,'sum')./exec_fun_on_cell_mat(allrespOR,'length');
    pallrespAND = 100*exec_fun_on_cell_mat(allrespAND,'sum')./exec_fun_on_cell_mat(allrespAND,'length');
    
    [mparOR,semparOR] = findMeanAndStandardError(pallrespOR);
    
    disp('Done');
    
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);

     %%


%% 1 off diagnoal (uniqe between adjacent trials)

for an = 1:5
    pcs.respRV(an,:) = diag(all_CI_mat(:,:,an));
    pcs.conjV(an,:) = diag(all_CI_mat(:,:,an),1);
    pcs.comp1V(an,:) = diag(uni(:,:,an),1);
    pcs.comp2V(an,:) = diag(uni(:,:,an),-1);
end
