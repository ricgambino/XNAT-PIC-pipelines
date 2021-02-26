#!/home/xnat/anaconda3/envs/py27/bin/python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 15:27:35 2019

@author: xnat
"""

# import matlab.engine

from REST_XNAT_Getscans_bytype import REST_XNAT_Getscans_bytype
from download_project_scans_of_type import download_project_scans_of_type
from process_scans import process_scans
from upload_XNAT_files import upload_XNAT_files
import sys

first_arg = sys.argv[1]
second_arg = sys.argv[2]

## default values are specified in order to test it without launching the pipeline on XNAT.
def process_T1_SR_project(
    build_dir="/data/xnat/build/",
    output_dir="/data/xnat/build/PROCESSED",
):
    #scan_type = ["epi","dwi"]
    scan_type=["rare-t1"]
    download_info = download_project_scans_of_type(build_dir, output_dir, scan_type)
    REST_XNAT_Getscans_bytype(download_info, build_dir)
    process_scans(download_info, build_dir)
    xnat_folder = "PROCESSED_T1_SR"
    upload_XNAT_files(download_info, build_dir, xnat_folder)


if __name__ == "__main__":
    
    process_T1_SR_project(first_arg, second_arg)
 
