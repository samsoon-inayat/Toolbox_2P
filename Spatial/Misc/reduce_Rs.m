function Rs = reduce_Rs(Rs,resp)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        Rs{rr,cc} = Rs{rr,cc}(resp{rr,cc},:);
    end
end

