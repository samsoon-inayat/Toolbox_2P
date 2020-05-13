function text = sReadTextFile (fileName)

fid = fopen(fileName,'r');
tline = fgetl(fid);
ii = 1;
text{ii} = tline;
while ischar(tline)
    text{ii} = tline;
    ii = ii + 1;
    tline = fgetl(fid);
end
fclose(fid);