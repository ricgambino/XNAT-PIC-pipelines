#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:21:06 2020

@author: Sara Zullino
"""

"""
Usage

mask_average('/path/to/T2/NifTi/file.nii', '/path/to/mask/NifTi/file.nii', '/output/path')

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

def mask_average(input_param_file , input_mask_file, output_dir):

    try:

        os.chdir(output_dir)

                
        print(input_mask_file)
        print(type(input_mask_file))
        print(input_mask_file.split("\n"))
        mask_files = input_mask_file.split("\n")
        # print("---")
        # print(mask_files[0])
        # print(mask_files[1])
        # print(len(mask_files))
        nROI = len(mask_files)
        
            
        # Load param image
        param_image = nib.load(input_param_file)
        size_x_param = param_image.shape[0]
        size_y_param = param_image.shape[1]
        try:
            size_z_param = param_image.shape[2]
        except:
            size_z_param = None    
        param_hdr = param_image.header
        param_qform = param_hdr['qform_code']
        param_sform = param_hdr['sform_code']
        canonical_param_image = nib.as_closest_canonical(param_image)
        param_img_orientation = nib.aff2axcodes(canonical_param_image.affine)
        param_image_data = param_image.get_fdata()
        param_image_original = param_image.get_fdata()
        param_filename = os.path.basename(input_param_file).split('.')[0]
        
        # Load mask
        mask_image = nib.load(mask_files[0])
        mask_image.affine.shape
        size_x_mask = mask_image.shape[0]
        size_y_mask = mask_image.shape[1]
        try:
            size_z_mask = mask_image.shape[2]
        except:
            size_z_mask = None 
        mask_hdr = mask_image.header
        mask_qform = mask_hdr['qform_code']
        mask_sform = mask_hdr['sform_code']
        canonical_mask = nib.as_closest_canonical(mask_image)
        mask_img_orientation = nib.aff2axcodes(canonical_mask.affine)
        mask = np.zeros((size_x_mask, size_y_mask, size_z_mask, nROI))
        for i in range(nROI):
            mask_image = nib.load(mask_files[i])
            mask_image_data = mask_image.get_fdata()
            mask[:,:,:,i] = mask_image_data
        
        
        # Check image orientation
        # print(param_hdr.get_sform(coded=True))
        # print(param_hdr.get_qform(coded=True))
        # print(mask_hdr.get_sform(coded=True))
        # print(mask_hdr.get_qform(coded=True))
        
        # mask_hdr.set_qform(param_hdr.get_sform, code='scanner')
        # param_hdr.set_qform(mask_hdr.get_sform, code='scanner')
        
        # if mask_qform == param_qform and mask_sform == param_sform:
        #     pass
        # elif mask_qform != param_qform and mask_sform != param_sform:
        #     mask_hdr.set_sform(param_hdr.get_sform, code='scanner')
        
        # affine = param_hdr.get_sform
        
        # img = nib.load(input_mask_file)
        # img.header.set_qform(affine, 1)
        # img.header.set_sform(affine, 1)
        # img.affine = affine
        # img_data = img.get_fdata()
        # plt.imshow(img_data,cmap='gray')
        # img.to_filename('tryout.nii')
        
        if size_x_mask != size_x_param and size_y_mask != size_y_param and size_z_mask != 1: 
            full_resized_mask = np.zeros((size_x_param, size_y_param, size_z_mask)) 
            with open("statistics.txt", "w+") as f:        
                f.write("***************************************************************************** \n")
                f.write("XNAT-PIC Pipeline: Mask Average\n") 
                f.write("Apply a mask to a parametric map and compute mean value in the ROIs\n")
                f.write("\n")
                f.write("Author: Sara Zullino \n")
                f.write("Mailto: sara.zullino@unito.it \n")    
                f.write("***************************************************************************** \n") 
                f.write("\n")   
                resized_mask = np.zeros((size_x_param,size_y_param,size_z_mask,nROI))
                for i in range(nROI):
                    mask_roi = mask[:,:,:,i]
                    #mask_copy_bin = np.sign(mask_copy)
                    mask_roi_resized = cv2.resize(mask_roi, dsize=(size_x_param, size_y_param), interpolation=cv2.INTER_LINEAR)    
                    mask_roi_resized_int = (mask_roi_resized != 0).astype(int)
                    resized_mask[:,:,:,i] = mask_roi_resized_int
                    temp = resized_mask[:,:,:,i]
                    full_resized_mask =  full_resized_mask + temp
              
                    # Applying multiply method 
                    param_image_mask = cv2.multiply(param_image_data, resized_mask[:,:,:,i], dtype=cv2.CV_32F)
                 
                    param_image_mask_nonzero_array = param_image_mask[param_image_mask>0]
                    mean_T2 = param_image_mask_nonzero_array.mean()
                    std_T2 = param_image_mask_nonzero_array.std()
                    median_T2 = np.median(param_image_mask_nonzero_array)  
                        
                    print("param_image_mask_nonzero_array max is" , np.amax(param_image_mask_nonzero_array))
                    print("param_image_mask_nonzero_array min is", np.amin(param_image_mask_nonzero_array))    
                
                    R2_image_mask_nonzero_array = 1000/param_image_mask_nonzero_array
                    mean_R2 = R2_image_mask_nonzero_array.mean()
                    std_R2 = R2_image_mask_nonzero_array.std()
                    median_R2 = np.median(R2_image_mask_nonzero_array)
                
                    # Count the number of nonzero pixels in the thresholded image
                    roi_area = cv2.countNonZero(resized_mask[:,:,:,i])
                    
                    print("ROI number: {}".format(i+1), file=f)
                    print("ROI file name: %s" % str(os.path.basename(mask_files[i])), file=f)
                    print("ROI Area: {:0.2f} ".format(roi_area), file=f)
                    print("Mean T2: {:0.2f} ms".format(mean_T2), file=f) 
                    print("STD T2: {:0.2f} ms".format(std_T2,2), file=f) 
                    print("Median T2: {:0.2f} ms".format(median_T2), file=f) 
                    print("Mean R2: {:0.2f} 1/s ".format(mean_R2,2), file=f) 
                    print("STD R2: {:0.2f} 1/s".format(std_R2,2), file=f) 
                    print("Median R2: {:0.2f} 1/s".format(median_R2), file=f) 
                    f.write("----------------------------------------------------------------------------- \n")        
            f.close()
                       
            # Save nifti figures  
            param_image_full_mask = np.zeros((size_x_param,size_y_param,size_z_mask))
            for j in range(0,size_z_mask):
                param_image_full_mask[:,:,j] = cv2.multiply(param_image_data[:,:,j], full_resized_mask[:,:,j], dtype=cv2.CV_32F)   
                #param_image_full_mask[:,:,j] = cv2.multiply(full_resized_mask[:,:,j], full_resized_mask[:,:,j], dtype=cv2.CV_32F) 
            full_resized_param_img = nib.Nifti1Image(param_image_full_mask, param_image.affine)
            full_resized_param_img.header.get_data_shape()
            nib.save(full_resized_param_img, '%s_masked.nii' % param_filename)
                 
            full_resized_mask_img = nib.Nifti1Image(full_resized_mask, mask_image.affine)
            full_resized_mask_img.header.get_data_shape()
            nib.save(full_resized_mask_img, 'mask.nii')    
            
            # # Figure: PARAMETRIC MAP
            # def display_multiple_img(images, rows = 1, cols=1):
            #     figure, ax = plt.subplots(nrows=rows,ncols=cols )
            #     for ind,title in enumerate(images):
            #         ax.ravel()[ind].imshow(images[title], cmap='gray')
            #         ax.ravel()[ind].set_title(title)
            #         ax.ravel()[ind].set_axis_off()
            #     plt.tight_layout()
            #     #plt.suptitle('%s masked [s]'% param_filename)
            #     plt.savefig('%s_masked.png' % param_filename)
            #     plt.show()
                 
            # total_images = size_z_mask
            # images = {str(i+1): param_image_full_mask[:,:,i] for i in range(total_images)}
            
            # display_multiple_img(images, 1, size_z_mask)
               
            # # Figure: MASK 
            # def display_multiple_img(images, rows = 1, cols=1):
            #     figure, ax = plt.subplots(nrows=rows,ncols=cols )
            #     for ind,title in enumerate(images):
            #         ax.ravel()[ind].imshow(images[title], cmap='gray')
            #         ax.ravel()[ind].set_title(title)
            #         ax.ravel()[ind].set_axis_off()
            #     plt.tight_layout()
            #     #plt.suptitle('mask')
            #     plt.savefig('mask.png')
            #     plt.show()
                 
            # total_images = size_z_mask
            # images = {str(i+1): full_resized_mask[:,:,i] for i in range(total_images)}
            
            # display_multiple_img(images, 1, size_z_mask)
             
        elif size_x_mask != size_x_param and size_y_mask != size_y_param and size_z_mask == 1:
            full_resized_mask = np.zeros((size_x_param, size_y_param))
            if size_x_mask != size_x_param and size_y_mask != size_y_param:
                resized_mask = np.zeros((size_x_param, size_y_param, nROI))
                with open("statistics.txt", "w+") as f:        
                    f.write("***************************************************************************** \n")
                    f.write("XNAT-PIC Pipeline: Mask Average\n") 
                    f.write("Apply a mask to a parametric map and compute mean value in the ROIs\n")
                    f.write("\n")
                    f.write("Author: Sara Zullino \n")
                    f.write("Mailto: sara.zullino@unito.it \n")    
                    f.write("***************************************************************************** \n") 
                    f.write("\n")
                    for i in range(nROI):
                        mask_roi = mask[:,:,:,i]
                        #mask_copy_bin = np.sign(mask_copy)
                        mask_roi_resized = cv2.resize(mask_roi, dsize=(size_x_param, size_y_param), interpolation=cv2.INTER_LINEAR)    
                        mask_roi_resized_int = (mask_roi_resized != 0).astype(int)
                        resized_mask[:,:,i] = mask_roi_resized_int  
                        temp = resized_mask[:,:,i]
                        full_resized_mask =  full_resized_mask + temp
        
                        # Applying multiply method 
                        param_image_mask = cv2.multiply(param_image_data, resized_mask[:,:,i], dtype=cv2.CV_32F)
                     
                        param_image_mask_nonzero_array = param_image_mask[param_image_mask>0]
                        mean_T2 = param_image_mask_nonzero_array.mean()
                        std_T2 = param_image_mask_nonzero_array.std()
                        median_T2 = np.median(param_image_mask_nonzero_array)              
                        
                        print("param_image_mask_nonzero_array max is" , np.amax(param_image_mask_nonzero_array))
                        print("param_image_mask_nonzero_array min is", np.amin(param_image_mask_nonzero_array))
                             
                        R2_image_mask_nonzero_array = 1000/param_image_mask_nonzero_array
                        mean_R2 = R2_image_mask_nonzero_array.mean()
                        std_R2 = R2_image_mask_nonzero_array.std()
                        median_R2 = np.median(R2_image_mask_nonzero_array)
                    
                        # Count the number of nonzero pixels in the thresholded image
                        roi_area = cv2.countNonZero(resized_mask[:,:,i])
                        
                        print("ROI number: {}".format(i+1), file=f)
                        print("ROI file name: %s" % str(os.path.basename(mask_files[i])), file=f)
                        print("ROI Area: {:0.2f} ".format(roi_area), file=f)
                        print("Mean T2: {:0.2f} ms".format(mean_T2), file=f) 
                        print("STD T2: {:0.2f} ms".format(std_T2,2), file=f) 
                        print("Median T2: {:0.2f} ms".format(median_T2), file=f) 
                        print("Mean R2: {:0.2f} 1/s ".format(mean_R2,2), file=f) 
                        print("STD R2: {:0.2f} 1/s".format(std_R2,2), file=f) 
                        print("Median R2: {:0.2f} 1/s".format(median_R2), file=f) 
                        f.write("----------------------------------------------------------------------------- \n")
                        
                f.close()
                
            # Save nifti figures  
            param_image_full_mask = cv2.multiply(param_image_data, full_resized_mask, dtype=cv2.CV_32F)   
         
            full_resized_param_img = nib.Nifti1Image(param_image_full_mask, param_image.affine)
            full_resized_param_img.header.get_data_shape()
            nib.save(full_resized_param_img, '%s_masked.nii' % param_filename)
                 
            full_resized_mask_img = nib.Nifti1Image(full_resized_mask, mask_image.affine)
            full_resized_mask_img.header.get_data_shape()
            nib.save(full_resized_mask_img, 'mask.nii')  
            
    except Exception as e:
        print(e)
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback)
        sys.exit(1)

if __name__=='__main__':
    mask_average(first_arg,second_arg,third_arg)    