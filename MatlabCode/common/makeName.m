function outName = makeName (name,fPath)
if iscellstr(name)
    outName = sprintf('%s\\%s',fPath,name{1});
else
    outName = sprintf('%s\\%s',fPath,name);
end