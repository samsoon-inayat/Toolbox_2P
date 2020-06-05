function T = check_abf_file(T)
for ii = 1:size(T,1)
    try
        rFolder  = cell2mat(T{ii,6});
    catch
        rFolder  = (T{ii,6});
    end
    files = dir(sprintf('%s\\*.abf',rFolder));
    if isempty(files)
        out(ii,1) = 0;
    else
        out(ii,1) = 1;
    end
end

TTemp = table(out,'VariableNames',{'ABF_File'});
T = [T TTemp];



