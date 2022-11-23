clear; clc; close all;

sizeCell = 50;
diameter = 20;
pxSize = 130;
thresh = 0.9;

dirPath = uigetdir;
fileList = dir([dirPath, filesep, '*_Simple Segmentation.h5']);
chns = ["_SNAP", "_RFP", "_mNG"]; chnsI = 1; % Wants to know distance between channel 1 and main channel

sol = [];
for n = 1 : length(fileList)
  % Loading in DAPI data and reshape to appopriate format
  rawData = h5read([dirPath, filesep, fileList(n).name], '/exported_data');
  rawData = reshape(rawData, 512, 512);
  rawData = rawData';

  % Get the channel used as identifier
  for m = 1 : length(chns)
    temp = strfind(fileList(n).name, chns(m));
    if ~isempty(temp)
      mainChns = m;
    end
  end

  % Calculating how many cells in the image that don't meet size requirement
  [cellData, cells] = bwlabel(rawData == 2);
  cells = [1:cells]';
  cellNumber = size(cells, 1);
  
  toRemove = [];
  for m = 1 : cellNumber
    if nnz(cellData == m) < sizeCell
      toRemove = [toRemove; m];
    end
  end
  cells(toRemove) = [];
  cellNumber = cellNumber - size(toRemove, 1);

  % Checking how many GFPs detected
  [gfpData, gfps] = bwlabel(rawData == 3);
  gfps = [1:gfps]';
  gfpNumber = size(gfps, 1);

  % Calculate distance between detected GFPs and SNAPs
  for m = 1 : cellNumber
    [x, y] = find(cellData == cells(m));
    xCell = x(boundary(x,y)); yCell = y(boundary(x,y));
    for i = 1 : gfpNumber
      [xGFP, yGFP] = find(gfpData == gfps(i));
      in = inpolygon(xGFP, yGFP, xCell, yCell);
      plot(y(boundary(x,y)),x(boundary(x,y)), 'LineWidth',2);
      hold on;
      scatter(yGFP, xGFP);
      set(gca, 'YDir', 'reverse')
      if sum(in) > 0
        imgs = zeros(diameter + 1, diameter + 1, length(chns));
        for j = 1 : length(chns)
          [stack, nbImages] = tiffread([dirPath, filesep, fileList(n).name(1:end-(23 + strlength(chns(mainChns)))), convertStringsToChars(chns(j)), '.tif']);

          % Convert to 3D matrix as double:
          imgs_3d_double = zeros(size(stack(1,1).data,1), size(stack(1,1).data,2), nbImages);
          for img_iter = 1:nbImages
            imgs_3d_double(:,:,img_iter) = stack(img_iter).data;
          end
          xGFP = floor((max(xGFP) - min(xGFP)) / 2) + min(xGFP);
          xGFP = floor(xGFP - (diameter / 2)) : ceil(xGFP + (diameter / 2));
          yGFP = floor((max(yGFP) - min(yGFP)) / 2) + min(yGFP);
          yGFP = floor(yGFP - (diameter / 2)) : ceil(yGFP + (diameter / 2));
          imgs(:, :, j) = imgs_3d_double(xGFP, yGFP);
          % Export the region of interest to analyse with different package if needed
%           imwrite(imgs(:,:,j)/ max(max(imgs(:,:,j))), [dirPath, filesep, fileList(n).name(1:end-(23 + strlength(chns(mainChns)))), convertStringsToChars(chns(j)), '_FOV.tif'])
        end
        imgI = imgs(:, :, chnsI) > (max(max(imgs(:, :, chnsI))) * thresh);
        imgI = imgs(:, :, chnsI) .* imgI;
        [imgI, dets] = bwlabel(imgI);
        [xT, yT] = find(imgs(:, :, mainChns) == max(max(imgs(:, :, mainChns))));
        xT = mean(xT); yT = mean(yT);
        mindist = (diameter + 1) * pxSize;
        for j = 1 : dets
          [xI, yI] = find(imgI == j);
          xI = mean(xI); yI = mean(yI);
          dist = sqrt((xT - xI) ^ 2 + (yT - yI) ^ 2) * pxSize;
          if dist < mindist
            mindist = dist;
          end
        end
        sol = [sol;  {fileList(n).name(1:end-(23 + strlength(chns(mainChns))))}, mindist];
%         imtool(imgs(:,:,1))
%         imshow(imgs(:,:,1) / max(max(imgs(:,:,1))))
      end
    end
  end
end
writecell(sol, [dirPath, filesep, 'Solution.csv']);