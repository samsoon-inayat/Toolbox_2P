try
    send_email({'samsoon.inayat@uleth.ca','brendan.mcallister@uleth.ca'},'Processing Protocol 10')
    ei10_C = loadContextsResponses(ei10_C,[1 1],[1 1 1]);
    send_email({'samsoon.inayat@uleth.ca','brendan.mcallister@uleth.ca'},'Done Protocol 10')
catch
    send_email({'samsoon.inayat@uleth.ca','brendan.mcallister@uleth.ca'},'Error occurred while loading data 10')
end

try
    send_email({'samsoon.inayat@uleth.ca','brendan.mcallister@uleth.ca'},'Processing Protocol 15')
    ei15_C = loadContextsResponses(ei15_C,[1 1],[1 1 1]);
    send_email({'samsoon.inayat@uleth.ca','brendan.mcallister@uleth.ca'},'Done Protocol 15')
catch
    send_email({'samsoon.inayat@uleth.ca','brendan.mcallister@uleth.ca'},'Error occurred while loading data 15')
end