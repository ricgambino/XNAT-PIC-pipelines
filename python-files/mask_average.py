#!/home/xnat/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  31 16:38:12 2020

@author: Sara Zullino
"""
"""
Usage

mask_average('/path/to/T2w/NifTi/file.nii', '/path/to/mask/NifTi/file.nii', '/output/path')

"""
import os
import sys
import traceback
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from statistics import median
import scipy 
import cv2

first_arg=sys.argv[1]
second_arg=sys.argv[2]
third_arg=sys.argv[3]

def mask_average(input_T2w_file , input_mask_file, output_dir):

    try:

        os.chdir(output_dir)

        # Load T2w image
        T2w_image = nib.load(input_T2w_file)
        size_x_T2w = T2w_image.shape[0]
        size_y_T2w = T2w_image.shape[1]
        T2w_hdr = T2w_image.header
        canonical_T2w_image = nib.as_closest_canonical(T2w_image)
        T2w_img_orientation = nib.aff2axcodes(canonical_T2w_image.affine)
        T2w_image_data = T2w_image.get_fdata()

        # Figure
        T2w_img_plot = plt.imshow(T2w_image_data,cmap='gray')
        plt.colorbar()
        plt.clim(0, 2)
        plt.title('T2w map [s]')
        plt.savefig('T2w_map.png')
        plt.clf()  
        plt.close()  

        # Load mask image
        mask_image = nib.load(input_mask_file)
        mask_image.affine.shape
        size_x_mask = mask_image.shape[0]
        size_y_mask = mask_image.shape[1]
        size_z_mask = mask_image.shape[2]
        mask_hdr = mask_image.header
        canonical_mask = nib.as_closest_canonical(mask_image)
        mask_img_orientation = nib.aff2axcodes(canonical_mask.affine)
        mask_image_data = mask_image.get_fdata()
        mask_image = mask_image_data[:,:,0]

        # Figure
        mask_img_plot = plt.imshow(mask_image,cmap='gray')
        plt.colorbar()
        plt.clim(0, 1)
        plt.title('mask') 
        plt.savefig('mask.png')
        plt.clf()  
        plt.close() 

        # Check image orientation
        #if mask_img_orientation != T2w_img_orientation:
            
        if size_x_mask != size_x_T2w and size_y_mask != size_y_T2w:
            mask_image_bin = scipy.sign(mask_image)
            mask_image_copy = mask_image.copy()
            mask_image_resized = cv2.resize(mask_image_bin, dsize=(size_x_T2w, size_y_T2w), interpolation=cv2.INTER_LINEAR)
            mask_image_resized_int = (mask_image_resized != 0).astype(int)
            plt.imshow(mask_image_resized_int,cmap='gray')
            
        # applying multiply method 
        T2w_image = np.asarray(T2w_image_data, dtype=None, order=None)
        T2w_image_masked = cv2.multiply(T2w_image, mask_image_resized_int,dtype=cv2.CV_32F)
        plt.imshow(T2w_image_masked,cmap='gray')
        plt.colorbar()
        plt.clim(0, 2)
        plt.title('T2w map masked [s]')
        plt.savefig('T2w_map_masked.png')
        plt.clf()  
        plt.close()

        # Count the number of nonzero pixels in the thresholded image
        roi_area = cv2.countNonZero(mask_image_resized_int)

        T2w_nonzero_array = T2w_image_masked[T2w_image_masked>0]
        mean_T2 = T2w_nonzero_array.mean()
        std_T2 = T2w_nonzero_array.std()
        median_T2 = np.median(T2w_nonzero_array)

        R2_nonzero_array = 1/T2w_nonzero_array
        mean_R2 = R2_nonzero_array.mean()
        std_R2 = R2_nonzero_array.std()
        median_R2 = np.median(R2_nonzero_array)


        with open("statistics.txt", "w+") as f:
            f.write("----------------------------------------------------------------------------- \n") 
            f.write("CIM-XNAT Pipeline: Tumor Masking \n") 
            f.write("Apply mask to a parametric map and compute mean value in the tumor region\n")
            f.write("Author: Sara Zullino \n") 
            f.write("----------------------------------------------------------------------------- \n") 
            print("ROI Area: {:0.2f} ".format(roi_area), file=f)
            print("Mean T2: {:0.2f} s".format(mean_T2), file=f) 
            print("STD T2: {:0.2f} s".format(std_T2,2), file=f) 
            print("Median T2: {:0.2f} s".format(median_T2), file=f) 
            print("Mean R2: {:0.2f} 1/s ".format(mean_R2,2), file=f) 
            print("STD R2: {:0.2f} 1/s".format(std_R2,2), file=f) 
            print("Median R2: {:0.2f} 1/s".format(median_R2), file=f) 
            f.close()

    except Exception as e:
        print(e)
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback)
        sys.exit(1)

if __name__=='__main__':
    mask_average(first_arg,second_arg,third_arg)    