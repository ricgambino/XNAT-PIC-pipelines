#!/home/xnat/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  28 16:36:12 2019

@author: Sara Zullino
"""
"""
Usage

from dcmrtstruct2nii import dcmrtstruct2nii, list_rt_structs

#print(list_rt_structs("/home/xnat/Documents/XNAT_pipelines_dev/apply_mask/aa_SL_rev_9_QG1_20181205_094817_untreated_post1w/AIM_20191105_171055/RTSTRUCT/AIM_20191105_171055.dcm"))

#dcmrtstruct2nii('/path/to/dicom/rtstruct/file.dcm', '/path/to/original/extracted/dicom/files', '/output/path')

"""
from dcmrtstruct2nii import dcmrtstruct2nii
import sys
import os
import traceback
import gzip, shutil

first_arg=sys.argv[1]
second_arg=sys.argv[2]
third_arg=sys.argv[3]

def rtstruct2nii(rtstruct_file, original_dcm, output_path):
    try:
        dcmrtstruct2nii(rtstruct_file,
                    original_dcm,
                    output_path)
    except Exception as e:
        print(e)
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback)
        sys.exit(1)
    try:
        with gzip.open(os.path.join(output_path,'mask_tumor.nii.gz'), 'r') as f_in, open(os.path.join(output_path,'mask_tumor.nii'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    except Exception as error:
        print(error)
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback)
        sys.exit(1)

if __name__=='__main__':
    rtstruct2nii(first_arg,second_arg,third_arg)    

		
		
		
		
		
		
		
