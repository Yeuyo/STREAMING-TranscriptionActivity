from pims import ImageSequence 
import trackpy as tp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

fileLoc = r"C:\Users\yeww\Downloads\STREAMING-Example\Image_data\ON\New folder (2)"
diameter = 5
cropSize = 20
pxSize = 130
locThresh = 390

rfpFiles = ImageSequence(fileLoc + "/*_RFP.tif")
mngFiles = ImageSequence(fileLoc + "/*_mNG.tif")
snapFiles = ImageSequence(fileLoc + "/*_SNAP.tif")
for n in range(len(rfpFiles)):
    # Get RFP image
    rfpImg = np.array(rfpFiles[n])
    rfpBG = rfpImg.mean()
    rfpSTD = rfpImg.std()
    # Get mNG image
    mngImg = np.array(mngFiles[n])
    mngBG = mngImg.mean()
    mngSTD = mngImg.std()
    # Get SNAP image
    snapImg = np.array(snapFiles[n])
    snapBG = snapImg.mean()
    snapSTD = snapImg.std()
    # Get RFP signals
    rfpSignals = tp.locate(rfpFiles[n], diameter)
    # Remove background noises
    rfpSignals = rfpSignals.loc[rfpSignals['signal'] > (2 * rfpBG),]
    rfpSignals.reset_index(drop = True, inplace = True)
    sol = pd.DataFrame()
    for m in range(len(rfpSignals)):
        # Calculate the location of the detected spot
        xCropL = math.floor(rfpSignals.loc[m, 'x'] - (cropSize / 2))
        xCropH = math.ceil(rfpSignals.loc[m, 'x'] + (cropSize / 2))
        yCropL = math.floor(rfpSignals.loc[m, 'y'] - (cropSize / 2))
        yCropH = math.ceil(rfpSignals.loc[m, 'y'] + (cropSize / 2))
        # Crop the spot out from all 3 channels
        rfpCrop = rfpImg[yCropL:yCropH, xCropL:xCropH]
        mngCrop = mngImg[yCropL:yCropH, xCropL:xCropH]
        snapCrop = snapImg[yCropL:yCropH, xCropL:xCropH]
        # Get the signals in the cropped image
        rfpCropSignal = tp.locate(rfpCrop, diameter, topn = 1)
        mngCropSignal = tp.locate(mngCrop, diameter, topn = 1)
        snapCropSignals = tp.locate(snapCrop, diameter)
        # Refining the signals
        rfpCropSignal = tp.refine_leastsq(rfpCropSignal, rfpCrop, diameter, fit_function='gauss')
        mngCropSignal = tp.refine_leastsq(mngCropSignal, mngCrop, diameter, fit_function='gauss')
        snapCropSignals = tp.refine_leastsq(snapCropSignals, snapCrop, diameter, fit_function='gauss')
        rfpCropSignal.reset_index(drop = True, inplace = True)
        mngCropSignal.reset_index(drop = True, inplace = True)
        # Show the detected spot
        plt.figure(0)
        plt.imshow(rfpCrop, vmin = rfpBG-rfpSTD*0.7, vmax = rfpBG+rfpSTD*4, cmap = 'gray')
        plt.scatter(rfpCropSignal['x'], rfpCropSignal['y'], marker="x", color = "blue", s = 250)
        plt.show()
        plt.figure(1)
        plt.imshow(mngCrop, vmin = mngBG-mngSTD*0.7, vmax = mngBG+mngSTD*4, cmap = 'gray')
        isLoc = False
        if mngCropSignal['signal'][0] > (2 * mngBG):
            plt.scatter(mngCropSignal['x'], mngCropSignal['y'], marker="x", color = "blue", s = 250)
            # Check for locality
            dist = math.sqrt((rfpCropSignal['x'][0] - mngCropSignal['x'][0])**2 + (rfpCropSignal['y'][0] - mngCropSignal['y'][0])**2)
            if (dist * pxSize) < locThresh:
                isLoc = True 
        plt.show()
        plt.figure(2)
        plt.imshow(snapCrop, vmin = snapBG-snapSTD*0.7, vmax = snapBG+snapSTD*4, cmap = 'gray')
        plt.scatter(snapCropSignals['x'], snapCropSignals['y'], marker="x", color = "blue", s = 250)
        plt.show()
        # Calculate minimum distance
        dist_min = []
        for Y, X in zip(snapCropSignals['y'], snapCropSignals['x']):
            distance = math.sqrt((Y - rfpCropSignal['y'][0])**2 + (X - rfpCropSignal['x'][0])**2)
            dist_min.append(distance)
        # sol.append([rfpFiles._filepaths[0].split('\\')[-1], isLoc, min(dist_min)* pxSize])
        sol = pd.concat([sol, pd.DataFrame({rfpFiles._filepaths[0].split('\\')[-1], isLoc, min(dist_min)* pxSize})], axis = 1)
sol = sol.T
sol.to_csv(fileLoc + "/Solution.csv")