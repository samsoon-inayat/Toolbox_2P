function conjunctive_representation

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0];
    rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets'};
    
    selContexts = [1 4 6 2 7 3 3 4 4 5 5 0 0];
    rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airIT','airD','airIT','airD','airIT','motionOnsets','motionOffsets'};
    
    o = get_data(ei,selContexts,rasterNames);
    for ii = 1:length(selContexts)
        xlabels{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    
    [speedRs,resp_speed] = get_speed_response(ei);
    xlabels{ii+1} = 'sp';
    resp = [o.resp.vals resp_speed];
    
    
    
    
    break
end
n = 0;
%%
an = 5;
[OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp);
maxOI = max([mOI(:);semOI(:)]); 
minOI = min([mOI(:);semOI(:)]);
figure(100);clf;
subplot 121
imagesc(mOI,[minOI maxOI]);
grid on;
set(gca,'xtick',1:length(xlabels),'ytick',1:length(xlabels),'xticklabels',xlabels,'yticklabels',xlabels,'Ydir','reverse'); xtickangle(45);colorbar;
axis equal
for rr = 1:size(OI_mat,1)
    for cc = 1:size(OI_mat,2)
        if isnan(OI_mat(rr,cc,1))
            continue;
        end
        if h_vals(rr,cc) == 1
            text(cc,rr,'*','Color','r');
        end
    end
end
subplot 122
imagesc(semOI,[minOI maxOI]);
grid on;
set(gca,'xtick',1:length(xlabels),'ytick',1:length(xlabels),'xticklabels',xlabels,'yticklabels',xlabels); xtickangle(45);colorbar;
axis equal


