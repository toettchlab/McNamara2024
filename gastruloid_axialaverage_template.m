% representative code analyze axial profiles of fluorescent markers in gastruloids
% N.B.: template modified and adapted for specific experiments
% written by H McNamara
% v1 April 2022
% updated Aug 2024


%% preliminaries
% make directories to store output
mkdir('normfigs');
mkdir('normdata');

% load in data 
close all;
clear all;
% takes in composite tif stack
fn = 'my_gastruloid_image';
imfn = [fn '.tif']
info = imfinfo(imfn);
numberOfPages = length(info)

oidnum = 0; %reset counter of segmented gastruloids
pagenumber = 0; % for loading multipate tiffs

%% read in and normalize image
% manually return to this block for each image within stack
offset = pagenumber*3;% to account for reading multipage tifs; adjust based on number of channels

iRFP = imread(imfn,1+offset);
GFP = imread(imfn,3+offset);
DsRed = imread(imfn,2+offset);
pagenumber = pagenumber + 1; %increment counter for next image

% visualize and confirm read in correcty
sz = size(DsRed);

sumImg = DsRed + GFP;
tmpRGB = zeros(sz(1), sz(2), 3);
tmpRGBscaled = zeros(sz(1), sz(2), 3);

tmpRGB(:,:,1) = DsRed;
tmpRGB(:,:,2) = GFP;
tmpRGB(:,:,3) = iRFP;

% quickly visualize with ad-hoc contrast adjustment
tmpRGBscaled(:,:,1) = double(DsRed)/1000;
tmpRGBscaled(:,:,2) = double(GFP)/1000;
tmpRGBscaled(:,:,3) = double(iRFP)/1000; 
imshow(tmpRGBscaled,[]);

% subtract background
% first, select area with no singals to get bkgd
[bkgdROI,bkgdvals] = clicky(tmpRGB, tmpRGBscaled); % using peripheral function which averages over user-selected area
% then subtract
DsRed2 = DsRed - bkgdvals(1);
GFP2 = GFP - bkgdvals(2);
iRFP2 = iRFP - bkgdvals(3);

% now, normalize images
tmp = sort(GFP2(:));
% take 99.9% percentile for max value
val = tmp(round(0.999*length(tmp))); 
% NB 1: for sparse labeling, increase threshold to avoid amplifcation weirdness
% NB 2: if no signal, can normalize based on hard-coded value from control issue
GFP3 = double(GFP2)/double(val);
imshow2(GFP3,[0 1]);

% repeat for other channels
tmp = sort(DsRed2(:));
val = tmp(round(0.999*length(tmp)));
DsRed3 = double(DsRed2)/double(val);
imshow2(DsRed3,[0 1]);

% % new in "_imagenorm" version: normalize image
tmp = sort(iRFP2(:));
val = tmp(round(0.999*length(tmp)));
iRFP3 = double(iRFP2)/double(val);
imshow2(iRFP3,[0 1]);

% display normalized image
tmpRGBscaled(:,:,1) = DsRed3;
tmpRGBscaled(:,:,2) = GFP3;
tmpRGBscaled(:,:,3) = iRFP3;
sumImg = GFP3 + DsRed3;
close all; imshow(tmpRGBscaled);

%% segment image
% first, binarize
imBin = imbinarize(sumImg);
imshow2(imBin,[]);
% fill in holes
imBin2 = imfill(imBin,'holes');
imshow2(imBin2,[]);

% dilate to fill out gaps near edges
seD = strel('disk',6,8);
imD = imdilate(imBin2,seD);
imshow2(imD,[]);

% erode away junk
seE = strel('diamond',8);
imE = imerode(imD,seE);
imshow2(imE,[]);

% optional step: manually seperate any touching gastruloids if necessary
close all;
[clipMask] = clicky_mask(imE);
tmp = imE - clipMask;
tmp2 = max(tmp,0);
imE = logical(tmp2);
% touch up to fill any residual holes
imF = imfill(imE,'holes');
imshow2(imF,[]);

% % extra manual step to fill holes if needed
% close all;
[fillMask] = clicky_mask(imF);
imF = or(imF, logical(fillMask));

% visualize to check segmentation quality
imRGB = zeros(sz(1), sz(2), 3);
imRGB(:,:,1) = sumImg*0.5;;
imRGB(:,:,3) = imF*0.5;
imshow(imRGB);

%% select gastruloid for analysis and generate A-P axial profiles
% select connected components for gastruloid of interest
close all; imshow2(imF,[]);
% note: for images with multiple gastruloids, manually repeat all of following code for each gastruloid
oidmask = bwselect(imF);
title('select oid for analysis')
oidnum = oidnum + 1; % increment counter

% check segementation
imRGB = zeros(size(tmpRGB));
imRGB(:,:,1) = DsRed3;
imRGB(:,:,2) = 0.5*oidmask;
imshow(imRGB);

%% skeletonize and average A-P profiles of signals

% skeletonize to find axis
close all;
oidskel = bwskel(oidmask,'MinBranchLength',200);
% note: MinBranchLength may require adjustment based on image
imshow(labeloverlay(10000*uint16(oidmask),oidskel,'Transparency',0))

% prune any spurious branches
close all;
tmpskelmask = zeros(size(tmpRGB));
tmpskelmask(:,:,1) = oidmask;
tmpskelmask(:,:,2) = oidskel;
[clipMask] = clicky_mask(tmpskelmask);
tmp = oidskel - clipMask;
tmp2 = max(tmp,0);
oidskel = logical(tmp2);
% check
imshow(labeloverlay(10000*uint16(oidmask),oidskel,'Transparency',0))

% manually extend skeleton to tips of gastruloids
[x, y] = meshgrid(1:sz(2), 1:sz(1));
% find location of endpoints
[i,j] = find(bwmorph(oidskel,'endpoints'));
i = [i(1); i(end)]; j = [j(1); j(end)]; % enforce only taking ends

for k = 1:2
    imshow(labeloverlay(10000*uint16(oidmask),oidskel,'Transparency',0));
    title('clicky to extend skeleton to ends')
    [xp, yp] = (getline(gca));
    enddists = sqrt((xp-j).^2 + (yp - i).^2);
    
    
    if enddists(1)<enddists(2) % find closest point to extend
        xv = [xp; j(1)];
        yv = [yp; i(1)];
    else
        xv = [xp; j(2)];
        yv = [yp; i(2)];
    end
    %     tmp = inpolygon(x,y,xv, yv);
        tmp = insertShape(uint16(oidskel),'Line',[xv yv]);
        oidskel = logical(tmp(:,:,1));
end

% re-skeleonize
oidskel2 = bwskel(oidskel);
% visualize to check
imshow(labeloverlay(10000*uint16(oidmask),oidskel2,'Transparency',0));

% specify posterior end by clickying close to point
% ideally, posterior is chosen based on a patterning readout (e.g. Wnt activity, Brachyury stain)
% without marker, can also specify based on morphology
close all;
imRGB2 = tmpRGBscaled .* repmat(oidmask,[1 1 3]);
imshow2(imRGB2,[]); title('right-click posterior end');
[xpost, ypost] = (getline(gca));
[i,j] = find(bwmorph(oidskel2,'endpoints'));
i = [i(1); i(end)]; j = [j(1); j(end)]; % enforce only taking ends
% further point is anterior end; closer is posterior
enddists = sqrt((xpost-j).^2 + (ypost - i).^2);

% compute distance along axis. 0 at anterior axis
pix = 1.4028; % pixel size in microns. adjust based on image resolution
if enddists(1)>enddists(2) %then first point is anterior pole
    D = bwdistgeodesic(oidskel2, j(1),i(1))*pix;
else
    D = bwdistgeodesic(oidskel2, j(end), i(end))*pix;
end
% display 
imshow2(D,[]);

% set up arrays for binning along axis
dL = 10; % linear bin size, in microns
oidL = max(D(:)); % overall axis length in microns
vecL = dL*(floor(oidL/dL)); %bins marked with starting value
L = 0:dL:vecL;
L_accum = zeros(size(L)); %set up count how many pixels assigned to each bin

% create axial coordinate mask of gastruloid
% i.e., mask which assigns each pixel within gastruloid to its closest
% position along A-P axis
[a,b] = find(oidskel2==1); % a,b points in row/column form
% now, find closest point for each point in segmented mask
oid_axmask = zeros(sz);
oid_binmask = zeros(sz);
for i = 1:sz(1)
    i
    for j = 1:sz(2) % for each pixel in image
        if oidmask(i,j) > 0 % if contained in segmentation
            % calculate distances
            dists_sqr = (i-a).^2 + (j-b).^2;
            % find closest point
            mindist = min(dists_sqr);
            idx = find(dists_sqr==mindist);
            % in case multiple points, take first in tiebreaker 
            idx = idx(1);
            % get distance along asks
            ap_coord = D(a(idx), b(idx));
            oid_axmask(i,j) = ap_coord;
            % find A-P profile bin to assign this coordinate to
            Lidx = floor(ap_coord/dL)+1;
            L_accum(Lidx) = L_accum(Lidx) + 1;
            oid_binmask(i,j) = Lidx;
        end
    end
end

% dislay
imshow2(oid_axmask,[])
imshow2(oid_binmask,[])

% now, use binmask average signals of interest
GFPprof = zeros(size(L));
DsRedprof = zeros(size(L));
iRFPprof = zeros(size(L));

for i = 1:sz(1)
    i
    for j = 1:sz(2) % for each pixel in image
        if oidmask(i,j) > 0 % if contained in segmentation
            Lidx = oid_binmask(i,j);
            % update signal accumulators, weighted by total counts at each
            % axial bin
            GFPprof(Lidx) = GFPprof(Lidx) + double(GFP3(i,j)) / L_accum(Lidx);
            DsRedprof(Lidx) = DsRedprof(Lidx) + double(DsRed3(i,j)) / L_accum(Lidx);
            iRFPprof(Lidx) = iRFPprof(Lidx) + double(iRFP3(i,j)) / L_accum(Lidx);
        end
    end
end

close all;

% generate figures
% first, normalized trace from all channels separately
close all; figure(); hold on;
plot(L,GFPprof,'g','LineWidth',1.5);
plot(L,DsRedprof,'r','LineWidth',1.5);
plot(L,iRFPprof,'b','LineWidth',1.5);
xlabel('distance along A-P axis (um)');
ylabel('normalized fluorescence (AU)');
legend('GFP','DsRed','iRFP','Location','best');
triplesave(gca, ['normfigs\' fn '_oid' num2str(oidnum) '_profiles']); % helper function which saves fig as .fig, .png, and .svg

% now, plot fraction of DsRed -> GFP labeling alongside iRFP
FracLabelProf = GFPprof ./ (GFPprof + DsRedprof);
close all; figure(); hold on;
plot(L,FracLabelProf,'k','LineWidth',1.5);
ylabel('GFP/(GFP+DsRed)')
legend('GFP/(GFP+DsRed)','location', 'best');
% plot(L,FracLabelProf);
triplesave(gca, ['normfigs\' fn '_oid' num2str(oidnum) '_fracprofile']);

% % now, save outputs
save(['normdata\' fn '_oid' num2str(oidnum) '_data.mat'], 'GFPprof', 'DsRedprof','iRFPprof', 'FracLabelProf','L','D','oidmask','oid_axmask', 'oid_binmask','oidskel','GFP','DsRed','iRFP', 'GFP3', 'DsRed3', 'iRFP3');

% repeat blocks as necessary to analyze all data


