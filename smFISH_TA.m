clear; clc; close all;

sizeCell = 50; % minimum size of cell to be considered as a cell
diameter = 20; % size of image to crop per detected mNG/RFP
pxSize = 130; % camera/miscrope pixel size
thresh = 0.9; % threshhold for computing blobs for SNAP
locDefine = 390; % RFP and mNG has to be less than this distance to be local

dirPath = uigetdir;
rfpList = dir([dirPath, filesep, '*_RFP.tif']);
mngList = dir([dirPath, filesep, '*_mNG.tif']);
snapList = dir([dirPath, filesep, '*_SNAP.tif']);

sol = [];
for n = 1 : length(rfpList)
  % Read RFP tif
  stack = tiffread([dirPath, filesep, rfpList(n).name]);
  imgR = stack(1).data;
  % Filter out background noise
  bgR = mean(mean(imgR));
  imgF = imgR > (10 * bgR);
  imgR = double(imgR) .* imgF;
  % Check for number of spots detected
  [imgI, dets] = bwlabel(imgR);

  % Read mNG tif
  stack = tiffread([dirPath, filesep, mngList(n).name]);
  imgM = stack(1).data;
  % Compute mNG background noise
  bgM = mean(mean(imgM));
  % Read SNAP tif
  stack = tiffread([dirPath, filesep, snapList(n).name]);
  imgS = stack(1).data;

  for m = 1 : dets
    [xdet, ydet] = find(imgI == m);
    % Determine the location of interest (place with RFPs)
    xCrop = floor((max(xdet) - min(xdet)) / 2) + min(xdet);
    xCrop = floor(xCrop - (diameter / 2)) : ceil(xCrop + (diameter / 2));
    yCrop = floor((max(ydet) - min(ydet)) / 2) + min(ydet);
    yCrop = floor(yCrop - (diameter / 2)) : ceil(yCrop + (diameter / 2));

    % Get center location of spot
    xdet = mean(xdet); ydet = mean(ydet);
    
    % Crop tifs
    imgRc = imgR(xCrop, yCrop);
    imgMc = imgM(xCrop, yCrop);
    imgSc = imgS(xCrop, yCrop);
    % Check for spots in RFP tif
    [xR, yR] = find(imgRc > 10 * bgR);
    xR = mean(xR); yR = mean(yR);
    % Check for spots in mNG tif
    [xM, yM] = find(imgMc > 10 * bgM);
    % Check for locality
    isLoc = 0;
    if ~isempty(xM)
      xM = mean(xM); yM = mean(yM);
      dist = sqrt((xR - xM) ^ 2 + (yR - yM) ^ 2) * pxSize;
      if dist < locDefine
        isLoc = 1;
      end
    end
    % Check for spots in SNAP
    imgSc = double(imgSc) .* (imgSc > (max(max(imgSc)) * thresh));
    [imgScI, detsSc] = bwlabel(imgSc);
    mindist = (diameter + 1) * pxSize;
    for i = 1 : detsSc
      [xS, yS] = find(imgScI == i);
      xS = mean(xS); yS = mean(yS);
      dist = sqrt((xR - xS) ^ 2 + (yR - yS) ^ 2) * pxSize;
      if dist < mindist
        mindist = dist;
      end
    end
    sol = [sol;  {rfpList(n).name(1:end-8)}, isLoc, mindist];
  end
end
writecell(sol, [dirPath, filesep, 'Solution.csv']);