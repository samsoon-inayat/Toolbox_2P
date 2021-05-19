function cn = get_context_number(ei,context_id_or_name)

CDs = contextDefinitions;
if ischar(context_id_or_name)
    selcName = context_id_or_name;
else
    selcName = CDs.names{context_id_or_name};
end

pp = 1;
contexts = ei.plane{1}.contexts;
for ii = 1:length(contexts)
    cNames{ii} = contexts(ii).name;
end
cn = find(strcmp(cNames,selcName));


