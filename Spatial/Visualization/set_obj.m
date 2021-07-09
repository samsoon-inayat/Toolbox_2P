function h = set_obj(ff,propsvals,vals)
if isstruct(ff)
    has = ff.h_axes;
else
    has = ff;
end

if ~exist('vals','var')
    for kk = 1:2:length(propsvals)
        set_prop(has,propsvals{kk},propsvals{kk+1});
    end
else
    set_prop(has,propsvals,vals);
end

function set_prop(has,prop,val)
if ~isequal(size(has),size(val))
    for ii = 1:size(has,1)
        for jj = 1:size(has,2)
            set(has(ii,jj),prop,val);
        end
    end
else
    for ii = 1:size(has,1)
        for jj = 1:size(has,2)
            set(has(ii,jj),prop,val{ii,jj});
        end
    end
end