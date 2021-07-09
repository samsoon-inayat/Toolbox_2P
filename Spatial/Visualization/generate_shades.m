function s = generate_shades(n)

shadg = linspace(0.4,1,n);
shad = linspace(0.4,1,n);

inds = (1:length(shad));
indsg = (1:length(shadg));

for ii = 1:length(shad)
    jj = inds(ii);
    s.m{jj,1} = [shad(ii) 0 shad(ii)];
    s.c{jj,1} = [0 shad(ii) shad(ii)];
    s.y{jj,1} = [shad(ii) shad(ii) 0];
end


for ii = 1:length(shadg)
    s.g{indsg(ii),1} = [shadg(ii) shadg(ii) shadg(ii)];
end



