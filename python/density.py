#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19, 2021

     This file is part of breast_density.

    breast_density is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    breast_density is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with breast_segmentation.  If not, see <http://www.gnu.org/licenses/>.

@author: Wolf-Dieter Vogl
"""
_LIMIT_BLOP_SIZE=1000  #minimum size of breast in #voxels
_BIN_SIZE=15 #2d histogram bin size to detect peaks

def _mergeFatWaterDecision(water_vol, fat_vol):
    """
    Merge fat-weighted and water-weighted image by computing 
    result=1./(1.+exp(-1.0*( water_vol-0.25)/fat_vol))

    water_vol and fat_vol are normalized to [0,1] before. 
    
    Parameters
    ----------
    water_vol : 3d volume
        volumetric information containing water-weighted values.
    fat_vol : TYPE
        volumetric information containing fat-weighted values.

    Returns
    -------
    result : 3d volume
        Merged volume

    """
    import numpy as np
    
    wn = np.float64(water_vol - np.min(water_vol))/np.float64(np.ptp(water_vol))
    fn = np.float64(fat_vol - np.min(fat_vol))/np.float64(np.ptp(fat_vol))
    
    result=1./(1.+np.exp(-1.0*( wn-0.25)/fn))
    
    return result

def _split_left_right(seg):
    """
    Split the breast segmentation in a left and right part. Splitting point is halfway between left and right breast bounding box determined
    by connected components. 
    
    If both breasts are connected, an alternative algorithm is used to determine splitting location. 
    Same for segmentation connected components containing less than _LIMIT_BLOP_SIZE pixel, assuming that these ones are only noise. 
    
    Parameters
    ----------
    seg : 3d binary volume
        Segmentation of both breasts.

    Returns
    -------
    seg_left, seg_right : 3d binary volumes.
        Segmentation of left and right breast only
    split_point : int
        x-coordinates halfway inbetween left and right segmentation
    """
    import numpy as np
    from skimage.measure import label, regionprops
    
    l = label(seg, connectivity=2)
    
    #compute regionprops and sort it by area
    r =  sorted(
            regionprops(l),
            key=lambda r: r.area,
            reverse=True,
    )
    
    #only one blop, use alternative
    if (len(r)<=1):
        return _segment_left_right(seg)
    
    #second blop is very small, probably noise only, use alternative
    if (r[0].area < _LIMIT_BLOP_SIZE or r[1].area < _LIMIT_BLOP_SIZE):
        return _segment_left_right(seg)
    
    #determine split point as middle distance between bounding box borders
    xmin = np.asarray([r[0].bbox[0],r[1].bbox[0]])
    xmax = np.asarray([r[0].bbox[3],r[1].bbox[3]])-1
    
    if (xmin[1] > xmin[0]):
        split_point = np.uint16(np.rint((xmin[1]-xmax[0])/2.0) + xmax[0])
    else:
        split_point = np.uint16(np.rint((xmin[0]-xmax[1])/2.0) + xmax[1])
    
    s_right = np.zeros_like(seg)
    s_left = np.zeros_like(s_right)
    s_right[:split_point+1,:,:] = seg[:split_point+1,:,:]
    s_left[split_point+1:,:,:] = seg[split_point+1:,:,:]
    
    return s_left,s_right,split_point

def _segment_left_right(seg):
    """Returns a segmentation of the breast into left and
    right breast when both segmetations are not separated.
    It computes a chessboard distance transformation from center slice of segmentation. 
    Distance is largest in the middle between the two breasts. 
    Maximum distance is detected by using minimum of second derivative
    
    Parameters
    ----------
    seg : 3d binary volume
        Segmentation of both breasts.

    Returns
    -------
    seg_left, seg_right : 3d binary volumes.
        Segmentation of left and right breast only
    split_point : int
        x-coordinates halfway inbetween left and right segmentation
    """
    
    import numpy as np
    from skimage.measure import label, regionprops
    from scipy.ndimage import distance_transform_cdt
    
    l = label(seg, connectivity=2)
    
    #compute regionprops and sort it by area, use only largest blop
    r =  sorted(
            regionprops(l),
            key=lambda r: r.area,
            reverse=True,
    )[0]
    
    xmin = r.bbox[0]
    ymin = r.bbox[1]
    zmin = r.bbox[2]
    
    xmax = r.bbox[3]
    ymax = r.bbox[4]
    zmax = r.bbox[5]
    
    d = np.float64(distance_transform_cdt(1-seg[xmin:xmax,ymin:ymax,zmin:zmax], metric='chessboard'))
    d = (d - np.min(d))/np.ptp(d)
    zmean = np.uint16(np.rint((zmax-zmin)/2.0))
    # Distance transformation has a 'w' like shape along y axis. The center of the
    #'w' shape is the center between the breasts.
    
    #Calculate 2nd derivatives to find "center" of w shape
    #Find minimum of 2nd derivative in the first 10 slices along y axis
    #Minimum is the x value between breasts
    #g = diff(d(:,1:10,zmean),2,1); 
    #Use central difference method to compute second derivative
    g = np.gradient(np.gradient(d[:,0:10,zmean],axis=0),axis=0)
    # find minimum in second derivative and x (column) coordinate
    x = np.unravel_index(np.argmin(g),np.shape(g))[0]
    split_point = x+xmin
    
    s_right = np.zeros_like(seg)
    s_left = np.zeros_like(s_right)
    s_right[:split_point+1,:,:] = seg[:split_point+1,:,:]
    s_left[split_point+1:,:,:] = seg[split_point+1:,:,:]
    return (s_left,s_right,split_point)


def _cluster_density_histogram (f, w, s):
    """
    Compute a 2d histogram of fat (f) weighted and water weighted (w) intensity pairs within segmentation (s). 
    Then a threshold is determined as halfway between histogram peaks. Finally, depending on the threshold for each voxel
    a label is assigned. 

    Parameters
    ----------
    f : 3d volume
        volume of fat weighted intensities
    w : TYPE
        volume of water weighted intensities
    s : TYPE
        segmentation of breast (breasts)

    Returns
    -------
    assignment : 3d volume
       label volume with background (0), parenchymal tissue (1), and fatty tissue (2) 
    histogram : 2d array
       2d histogram of fat and water intensities
    histogram_threshhold: 2d array
       2d label image with same size than histogram indicating the threshold areas. 
    peak-center : 2 element list
        center position of histogram peaks
    
    """
    import numpy as np
    from skimage.morphology import h_maxima,label
    from skimage.measure import regionprops
    
    #get intensity values within segmentation only and create pairs
    X = np.array([f[s>0],w[s>0]],dtype=np.uint32)
    #get maximum intensities
    maxf = np.uint32(np.max(X[0,:]))+1
    maxw = np.uint32(np.max(X[1,:]))+1
    
    #Calculate histogram and binned histogram 
    #A=np.zeros((maxf,maxw),dtype=np.uint32)
    #for c in X.T: #iterate over each pair and count up in histogram matrix
    #    A[c[0],c[1]]+=1
    A,_,_ = np.histogram2d(X[0,:],X[1,:],bins=[np.arange(maxf+1),np.arange(maxw+1)])
    A = np.uint32(A)
    N,_,_ = np.histogram2d(X[0,:],X[1,:],bins=_BIN_SIZE)
    
    #Find extrema in binned histogram    
    J = h_maxima(N,500)
    L = label(J)
    numL = np.max(L)
    stats = regionprops(L)
    
    value = []
    #get sorting of peaks by histogram count and filter invalid peaks 
    for stat in stats:
        #often a peak is at 0,0 when segmentation is outside of breast area, where both intensities are close to 0
        #ignore this peak
        if (np.max(stat.centroid)==0):
            value.append(-np.inf)
        else:
            centroid = np.uint16(np.rint(stat.centroid))
            value.append(N[centroid[1],centroid[0]])
        #sort by highest value
        idx = np.flip(np.argsort(value))
            
    #first peak center is the one with highest value
    bincenter1 = np.array(stats[idx[0]].centroid)
    
    #initialize second peak center with mirrored position along diagonal
    bincenter2 = np.array([bincenter1[1],bincenter1[0]])
    
    
    if (numL > 2):
        #more than two peaks, choose peaks with highest histogram values and being on opposite site
           
        #iterate over all peaks and choose the one that is on opposite site
        for i in idx[1:]:
            if (value[i] == -np.inf):
                continue
            t = stats[idx[i]].centroid
            #check if it is on opposite site
            if ((bincenter1[0] >= bincenter1[1] and t[0] <= t[1]) or
                (bincenter1[0] <= bincenter1[1] and t[0] >= t[1])):
                bincenter2 = t
                break

    elif (numL==2):
        #Check whether detected peaks are on opposite site
        if ((bincenter1[0] >= bincenter1[1] and bincenter2[0] <= bincenter2[1]) or
            (bincenter1[0] <= bincenter1[1] and bincenter2[0] >= bincenter2[1])):
            # cluster centers are not on different sides, peak with highest value is chosen and mirrored
                bincenter2 = stats[1].centroid
    
    #get interval size
    h = np.float64(np.max(X,axis=1)) / _BIN_SIZE
    #convert coordinates from binned histogram [0..BIN_SIZE) to original histogram coordinates (center of interval = (h/2) is used)
    center1 = bincenter1 * h + h/2.0
    center2 = bincenter2 * h + h/2.0
    
    #determine halfway position between the two peak centers
    center = (center1 + center2)/2.0
    
    # create splitting matrix depending on center position
    a = center[1]/center[0]
    s1 = np.int32(a * np.tile(np.arange(A.shape[0]),[A.shape[1],1]).T)
    s2 = np.tile(np.arange(A.shape[1]),[A.shape[0],1])
    s_split = s1-s2
    s_split = s_split > 0
    
    #s_split has now same size as 2d-histogram and contains information whether each element is fatty tissue or water tissue
    
    #now perform a lookup for each voxel within segmentation whether it is fatty or water tissue    
    #lookup indices are the paired intensity values in s_split
    assignment = 1 + s_split[tuple(X)] #  + 1 to get label 1/2 (fat/water)
    #reshape assignment to segmented breast
    assign_vol = np.zeros_like(s,dtype=np.uint8)
    assign_vol[s>0] = assignment
    
    return (assign_vol, A, s_split, [center1,center2])
    
  
def _compute_density_values (assignment, voxel_spacing,postfix=''):
    """
    Compute from a fat/water segmentation some summary measures (total volume, fat volume, water volume, and %)

    Parameters
    ----------
    assignment : 3d label volume
        volume containing laber information for background, fat and water tissue.
    voxel_spacing : 3-element vector
        Spacing information of voxel in mm for all 3 dimensions
    postfix : string
        this string is added as postfix to each dictionary key. This is usefull to get distinguish for instance between left and right breast
        in the results
        
    Returns
    -------
    density : dictionary
        volume              Volume of both breasts [cm³]
        volume_fat          Volume of fatty tissue in both breasts [cm³]
        volume_water        Volume of parenchymal tissue in both breasts [cm³]
        percentage_fat           Percentage of fatty tissue in breasts [%]
        percentage_water            Percentage of parenchymal tissue in breasts [%]
    """    
    import numpy as np
    density = {}
    density['volume_fat'+postfix]= np.sum(assignment==1) * voxel_spacing[0] * voxel_spacing[1] * voxel_spacing[2] * 1000.0 #convert from pixel spacing to physical spacing in mm³ to cm³
    density['volume_water'+postfix]= np.sum(assignment==2) * voxel_spacing[0] * voxel_spacing[1] * voxel_spacing[2] * 1000.0 
    density['volume'+postfix]= density['volume_fat'+postfix] + density['volume_water'+postfix]
    density['percentage_fat'+postfix]= density['volume_fat'+postfix] / density['volume'+postfix] * 100  #%
    density['percentage_water'+postfix]= density['volume_water'+postfix] / density['volume'+postfix] * 100 #%
    
    return density

def _plot_results(w, assign_vol, split_point, hist, hist_split, center_pos, output_folder):    
    """
    Plot fat/water assignment, split left/right breast, histogram and histogram split
     
    Parameters
    ----------   
    w : 3d volume
        water weighted volume, that is used for overlaying segmentations
    assign_vol : 3d label volume
        label for fat/water assignment
    split_point : scalar
        x-position of left-right breast separation
    hist : 2d array
        2d histogram that was used to determine optimal threshold
    hist_split : 2d array
        threshold for histogram
    center_pos : 2 element list
        center position of histogram peaks
    output_folder: string or Path
        location where figures should be saved
    """
    
    #get bounding box and axial center slice
    from skimage.measure import regionprops
    import matplotlib.pyplot as plt
    from pathlib import Path
    import numpy as np
    
    segname = Path(output_folder).joinpath('div_left_right.png')
    assignname = Path(output_folder).joinpath('assignment.png')
    histname = Path(output_folder).joinpath('hist.png')
    histdivname = Path(output_folder).joinpath('hist_div.png')
    
    s = np.ones_like(assign_vol,dtype=np.float32)
    s[assign_vol==0]=np.NaN
    
    bbox = regionprops(np.uint8(assign_vol>0))[0].bbox
    xmin = bbox[0]
    xmax = bbox[3]
    ymin = bbox[1]
    ymax = bbox[4]
    zmean = np.uint16(np.rint((bbox[5]-bbox[2])/2.0))
    
    #plot segmentation
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.imshow (w[xmin:xmax,ymin:ymax,zmean].T,cmap=plt.cm.gray,interpolation='nearest')
    ax.imshow (s[xmin:xmax,ymin:ymax,zmean].T,cmap=plt.cm.Reds,alpha=.75,vmin=0,vmax=1, interpolation='nearest')
    ax.plot([split_point-xmin,split_point-xmin],[0,ymax-ymin],'r-')
    ax.axis('off')
    fig.savefig(str(segname))
    plt.close(fig)
    
    #plot assignment
    assign = np.zeros_like(assign_vol,dtype=np.float32)
    assign[assign_vol==0]=np.NaN
    assign[assign_vol==1]=2 #revert to better match color map
    assign[assign_vol==2]=1
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.imshow (w[xmin:xmax,ymin:ymax,zmean].T,cmap=plt.cm.gray,interpolation='nearest')
    ax.imshow (assign[xmin:xmax,ymin:ymax,zmean].T,cmap=plt.cm.tab10,alpha=.75,vmin=0,vmax=2,interpolation='nearest')
    ax.axis('off')
    fig.savefig(str(assignname))
    plt.close(fig)   
    
    #plot histogram
    fig = plt.figure()
    ax = fig.add_subplot()
    
    ax.imshow (np.log(hist.T+1e-6),cmap=plt.cm.magma,origin='lower',aspect='equal')
    ax.set_ylabel('water')
    ax.set_xlabel('fat')     
    fig.savefig(str(histname))
    plt.close(fig)     
    
    #plot histogram assignment
    fig = plt.figure()
    ax = fig.add_subplot()

    ax.imshow (np.log(hist.T+1e-6),cmap=plt.cm.gray,origin='lower',aspect='equal')
    ax.imshow (hist_split.T,cmap=plt.cm.Set1,origin='lower',alpha=.5,vmin=0,vmax=1, interpolation='nearest',aspect='equal')
    ax.set_ylabel('water')
    ax.set_xlabel('fat')  
    fig.savefig(str(histdivname))
    plt.close(fig)       
    
    
def measure_density (water_img_filename, fat_img_filename, breast_segmentation_filename,  output_folder=None, save_plots=False,):
    """
    Compute a segmentation of parenchym and fatty tissue based on MRI images containing water weighted / fat weighted information, as
    the dixon sequence provides it, or alternatively any fat-suppressed/non-fat-suppressed weighting images. 
    All images are assumed to be in RL / AP / IS orientation and are already co-registered and having same pixel spacing. 
    Intensity values are processed as int16. 
    Orientation can be changed with fsl tool:
    fslswapdim filein.nii.gz RL AP IS fileout.nii.gz

    Parameters
    ----------
    water_img_filename : string or Path
        Path to nifti image file containing water weighted information
    fat_img_filename : string or Path
        Path to nifti image file containing fat weighted information
    breast_segmentation_filename : string or Path
        Path to nifti image file containing segmentation of breast (1 breast, all other values background)
    output_folder : string or Path
        optional. Folder where results are saved. Results are "density.csv" containing measurements about density and "assignment.nii" containing
        the labels for background (0), parenchymal tissue (1), and fatty tissue (2) 
        If None, no results are saved. 
        Default: None
    save_plots : boolean
        optional. When true, plots of histogram (hist.png), histogram threshold (hist_div.png), central axial slice
        illustrating division into left/right breast (div_left_right.png), and central axial slice of assignment (assignment.png)
        Parameter "output_folder" needs to be set, as in that location the results are saved. 
        Default: False
    Returns
    -------
    density : dictionary
        volume              Volume of both breasts [cm³]
        volume_left         Volume of left breast [cm³]
        volume_right        Volume of right breast [cm³]
        volume_fat          Volume of fatty tissue in both breasts [cm³]
        volume_fat_left/right     Volume of fatty tissue in left/right breast [cm³]
        volume_water        Volume of parenchymal tissue in both breasts [cm³]
        volume_water_left/right        Volume of parenchymal tissue in left/right breast [cm³]
        percentage_fat           Percentage of fatty tissue in breasts [%]
        percentage_fat_left/right   Percentage of fat in left/right breast [%]
        percentage_water            Percentage of parenchymal tissue in breasts [%]
        percentage_water_left/right Percentage of Water in left/right breast [%]
        split_point      Splitting point (x pixel coords) of left and right breast    
    assignment : nibabel nifti object
        volumetric segmentation as nifti object containing labels for background (0), parenchymal tissue (1), and fatty tissue (2)        
    """
    import nibabel as nib    
    from pathlib import Path
    import numpy as np
    from scipy.ndimage import binary_fill_holes
    import pandas as pd
    
    #Parameter verification
    #----------------------
    if save_plots and output_folder is None:
        raise ValueError("output_folder must be set when save_plots=True")
    

        
    #loading and preprocessing
    #-------------------------
    water_nii = nib.load(str(water_img_filename))
    fat_nii = nib.load(str(fat_img_filename))
    seg_nii = nib.load(str(breast_segmentation_filename))
    
    # we work with 16 bit integer data
    # --------------------------------
    w = np.int32(water_nii.get_fdata())
    f = np.int32(fat_nii.get_fdata())

    #combine w+f and normalize it between [0,1]
    #-----------
    fw = np.float64(w+f)        
    fw = (fw - np.min(fw))/np.ptp(fw)
    
    #combine w and f using a decision function
    wfdec = _mergeFatWaterDecision(w,f)
    
    #correct breast segmentation errors (segmentation of air-like structure)
    #-------------------
    # everything above 0.5 and below 1.5 is assumed as breast label
    s = np.logical_and(seg_nii.get_fdata()>=0.5,seg_nii.get_fdata()<1.5)
    #remove errornous segmentations of airlike structures (empirical wfdec < 0.01 or f+w < 0.05)
    s [wfdec<=0.01]=0
    s [fw<=0.05]=0
    s = binary_fill_holes(s)
    
    #split into left / right
    s_left, s_right, split_point = _split_left_right(s)

    #compute segmentation
    assign_vol, hist, hist_split, center_pos = _cluster_density_histogram (f,w,s)
    assign_left = assign_vol.copy()
    assign_right = assign_vol.copy()
    assign_left[s_left==0]=0
    assign_right[s_right==0]=0

    #compute density 
   
    density = _compute_density_values(assign_vol,seg_nii.header.get_zooms())
    density.update(_compute_density_values(assign_left,seg_nii.header.get_zooms(),'_left'))
    density.update(_compute_density_values(assign_right,seg_nii.header.get_zooms(),'_right'))
    
    #plot if activated
    if (save_plots):
        output_path = Path(output_folder)
        output_path.mkdir(parents=True,exist_ok=True)
        _plot_results(w, assign_vol, split_point, hist, hist_split, center_pos, output_folder)
    
    assign_nii = nib.Nifti1Image(assign_vol, water_nii.affine)
    #save results and plots
    if (output_folder is not None):
        output_path = Path(output_folder)
        output_path.mkdir(parents=True,exist_ok=True)
        df = pd.DataFrame(data=density,index=[1])
        df.to_csv(output_path.joinpath('density.csv'),index=False)
        nib.save(assign_nii,output_path.joinpath('assignment.nii'))
    
    #return
    return (density,assign_nii)

if __name__ == "__main__":
    measure_density('/home/wvogl/temp/pinker/w_t1.nii.gz','/home/wvogl/temp/pinker/f_t1.nii.gz','/home/wvogl/temp/pinker/s_wfdec.nii.gz', '/home/wvogl/temp/pinker/output',True)
    
