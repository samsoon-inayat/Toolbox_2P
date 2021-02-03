function out = get_mfile_vars(fileName)

run(fileName);
listofvars = who;

for ii = 1:length(listofvars)
    if ~isempty(strfind(listofvars{ii},'context'))
        cmdTxt = sprintf('out.%s = %s;',listofvars{ii},listofvars{ii});
        eval(cmdTxt);
    end
end
