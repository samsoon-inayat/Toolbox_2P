%%
x = triu(ones(size(aCRc{1,1})),1).*aCRc{1,1};
vals = x(mask);
figure(2000);clf;hist(vals,2000);
% x = zeros(size(aCRc{1,1}));
% x(:,60:76) = 1;
% x(60:76,:) = 1;
figure(1000);clf;
subplot 131
imagesc(x);axis equal; axis off;
colorbar;
xdft=fft2(x);
oxdft = xdft;
xdft=xdft.*conj(xdft);
xdft=xdft/(size(x,1)*size(x,2));
xdft=fftshift(xdft);
subplot 132
imagesc(xdft);
axis equal; axis off
subplot 133
imagesc(angle(oxdft));
axis equal; axis off

%%
tmR = mR{1,1};

[pc.coeff,pc.score,pc.latent,pc.tsquared,pc.explained,pc.mu] = pca(tmR);
figure(3000);clf;
imagesc(pc.score);