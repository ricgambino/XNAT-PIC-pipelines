#!/home/xnat/anaconda3/envs/py27/bin/python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 16:32:35 2019

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
def process_EPI_project(
    build_dir="/data/xnat/build/MITIGATE/20190620_222403/lc_41_IMA_T1_UNTREATED_PRE/RAW",
    output_dir="/data/xnat/build/MITIGATE/20190620_222403/lc_41_IMA_T1_UNTREATED_PRE/PROCESSED",
):
    scan_type=["dwi","epi"]
    download_info = download_project_scans_of_type(build_dir, output_dir, scan_type)
    REST_XNAT_Getscans_bytype(download_info, build_dir)
    process_scans(download_info, build_dir)
    xnat_folder = "PROCESSED_DWI"
    upload_XNAT_files(download_info, build_dir, xnat_folder)


if __name__ == "__main__":
    
    process_EPI_project(first_arg, second_arg)
    #process_EPI_project() ## To launch debugging
