function mask = getMask (dataFolder,varargin)

files = dir(dataFolder);

mi = [];
for ii = 1:length(files)
    if ~files(ii).isdir
        pos = strfind(lower(files(ii).name),'mask');
        if ~isempty(pos)
            mi = [mi ii];
        end
    end
end

if isempty(mi)
    mask = [];
    return;
end

for ii = 1:length(mi)
    pos = strfind(files(mi(ii)).name,'tif');
    if ~isempty(pos)
        break;
    end
end

fileName  = makeName(files(mi(ii)).name,dataFolder);
mask = tif2mat(fileName);
if max(mask) > 1
    umask = unique(mask(:));
    if length(umask) > 2
        error;
    else
        mask = mask==umask(1);
    end
end

mask = logical(mask);

if nargin == 2
    ec = varargin{1};
else
    return;
end

temp = bwboundaries(mask);
mask = false(size(mask));
for p = 1:length(temp)
    boundary = fliplr(temp{p});
    [ geom, iner, cpmo ] = polygeom( boundary(:,1), boundary(:,2) );
    x_cent = geom(2); y_cent = geom(3);
    cent = x_cent + 1i * y_cent;
    for ii = 1:size(boundary,1)
        thisP = boundary(ii,1) + 1i * boundary(ii,2);
        diffV = thisP - cent;
        l_diffV = abs(diffV);
        nl_diffV = ec * l_diffV;
        ang = angle(diffV);
        nThisP = (nl_diffV * cos(ang) + 1i *nl_diffV*sin(ang)) + cent;
        nBoundary(ii,1) = real(nThisP); nBoundary(ii,2) = imag(nThisP);
    end
    mask = mask + poly2mask(nBoundary(:,1),nBoundary(:,2),size(mask,1),size(mask,2));
end
mask = logical(mask);




