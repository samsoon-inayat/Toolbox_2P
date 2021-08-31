function print_names(Rs)

for rr = 1%:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        disp(R.marker_name);
    end
end