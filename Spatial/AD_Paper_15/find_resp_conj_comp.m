function [respRV,conjV,comp1V,comp2V] = find_resp_conj_comp(all_resp)

[OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_resp,0.5,0.05);

for an = 1:size(all_resp,1)
    respRV(an,:) = diag(all_CI_mat(:,:,an));
    conjV(an,:) = diag(all_CI_mat(:,:,an),1);
    comp1V(an,:) = diag(uni(:,:,an),1);
    comp2V(an,:) = diag(uni(:,:,an),-1);
end