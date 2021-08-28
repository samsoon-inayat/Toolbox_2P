function conjunctive_representation

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0];
    rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets'};
    o = get_data(ei,selContexts,rasterNames);
    
    selContexts = [3 4 5];
    rasterNames = {'air77T','air77T','air77T'};
    opct = get_data(ei,selContexts,rasterNames);
    
    for ii = 1:length(selContexts)
        xlabels{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    
    [speedRs,resp_speed] = get_speed_response(ei);
    xlabels{ii+1} = 'sp';
    resp = [o.resp.vals resp_speed];
%     resp_o = cell_list_op(resp,[],'xor')
    si_light = [1 2 3];
    si_air_rest = [4 5];
    si_air_dist_trials = [6 8 10];
    si_air_dist_itrials = [7 9 11];
    si_air_time_trials = [12 14 16];
    si_air_time_itrials = [13 15 17];
    si_motion = [18 19];
    
    dzMI = prop_op(o.props.zMI(:,[si_air_dist_trials si_air_dist_itrials]),o.props.zMI(:,[si_air_time_trials si_air_time_itrials]),0.1);
    
    resp = [resp dzMI.resp_D_g_T(:,1) dzMI.resp_T_g_D(:,1)];
    resp_OR = cell_list_op(resp,[],'or');
    resp_AND = cell_list_op(resp,[],'and');
    xlabels = [xlabels {'D_g_T','T_g_D'}];
    
    break
end
n = 0;
%%
si = si_air_time_trials;
all_conds = 1:size(resp,2);
resp_o = get_cell_list(resp,[si]');
resp_o1 = get_cell_list(resp,[-setdiff(all_conds,si)]);
resp_ = cell_list_op(resp_o,resp_o1,'and');
view_population_vector(o.Rs(:,si),o.mR(:,si),resp_o(:,si),100);
%%
an = 5;
[OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
maxOI = max([mOI(:);semOI(:)]); 
minOI = min([mOI(:);semOI(:)]);
maxOI = max(OI_mat,[],3);
figure(100);clf;
% subplot 121
fOI = maxOI;
imagesc(fOI,[min(fOI(:)) max(fOI(:))]);
grid on;
set(gca,'xtick',1:length(xlabels),'ytick',1:length(xlabels),'xticklabels',xlabels,'yticklabels',xlabels,'Ydir','reverse'); xtickangle(45);colorbar;
axis equal
for rr = 1:size(OI_mat,1)
    for cc = 1:size(OI_mat,2)
        if isnan(OI_mat(rr,cc,1))
            continue;
        end
        if h_vals(rr,cc) == 1
            text(cc,rr,'*','Color','r','FontSize',12);
        end
    end
end
% subplot 122
% imagesc(semOI,[minOI maxOI]);
% grid on;
% set(gca,'xtick',1:length(xlabels),'ytick',1:length(xlabels),'xticklabels',xlabels,'yticklabels',xlabels); xtickangle(45);colorbar;
% axis equal


