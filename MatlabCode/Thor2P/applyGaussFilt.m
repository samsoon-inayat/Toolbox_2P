function filtered = applyGaussFilt(sig,win)
GF = gausswin(win);
GF = GF/sum(GF);
dSigT = [zeros(1,win) sig];
dSigT = filter(GF,1,dSigT);
filtered = dSigT(win+1:end);