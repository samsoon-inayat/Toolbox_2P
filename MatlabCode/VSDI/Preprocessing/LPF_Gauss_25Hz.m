function H = LPF_Gauss_25Hz
H = gaussmf(-4:4,[1.1 0]) ;
H = H/sum(H);