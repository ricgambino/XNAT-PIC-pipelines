#!/home/xnat/anaconda3/envs/py27/bin/python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 16:29:19 2019

@author: xnat
"""
import os

# import subprocess
import matlab.engine


def process_scans(download_info, path):
    eng = matlab.engine.start_matlab()
    for j in range(0, len(download_info)):
        print(download_info)
        for scan_dirs in download_info[j]["Scan Dirs"]:           
            if not os.path.exists(download_info[j]["Subject Folder"]):
                os.mkdir(download_info[j]["Subject Folder"])
            if not os.path.exists(
                download_info[j]["Subject Folder"]
                + "/"
                + download_info[j]["Session ID"]
            ):
                os.mkdir(
                    download_info[j]["Subject Folder"]
                    + "/"
                    + download_info[j]["Session ID"]
                )
            if not os.path.exists(
                download_info[j]["Subject Folder"]
                + "/"
                + download_info[j]["Session ID"]
                + "/"
                + download_info[j]["Scan ID"]
            ):
                os.mkdir(
                    download_info[j]["Subject Folder"]
                    + "/"
                    + download_info[j]["Session ID"]
                    + "/"
                    + download_info[j]["Scan ID"]
                )
            try:
                # pass
                # print("Sono in process scan, folder is "+scan_dirs)
                # print("output folder is "+download_info[j]["Subject Folder"]+'/'+download_info[j]['Session ID'])
                # print(download_info[j]["Subject Folder"]+'/'+download_info[j]['Session ID'])
                eng.cd("/home/xnat/Documents/XNAT-PIC/XNAT-PIC-pipelines/stable/matlab-files", nargout=0)
                # print(eng.ls())
	        print(download_info[j]["Scan Type"])
		if download_info[j]["Scan Type"] == ['EPI diffusion map']:
                    eng.process_DWI(
                    scan_dirs,
                    download_info[j]["Subject Folder"]
                    + "/"
                    + download_info[j]["Session ID"]
                    + "/"
                    + download_info[j]["Scan ID"],
                    nargout=0,
                )
		elif download_info[j]["Scan Type"] == ['MSME-T2-map in vivo short'] or download_info[j]["Scan Type"] == ['MSME-T2-map in vivo']:
                    eng.process_T2w(
                    scan_dirs,
                    download_info[j]["Subject Folder"]
                    + "/"
                    + download_info[j]["Session ID"]
                    + "/"
                    + download_info[j]["Scan ID"],
                    nargout=0,
                )
                elif download_info[j]["Scan Type"] == ['RARE-T1 map vivo MTX64']:
                    eng.process_T1_SR(
                    scan_dirs,
                    download_info[j]["Subject Folder"]
                    + "/"
                    + download_info[j]["Session ID"]
                    + "/"
                    + download_info[j]["Scan ID"],
                    nargout=0,
                )

            except Exception as err:
                print(err)
                #eng.quit()
            finally:
                eng.quit()
    #eng.quit()
    # subprocess.check_call([script_path, src+scan_dirs, download_info[j]["Subject Folder"]+'/'+scan_dirs+'/'])
