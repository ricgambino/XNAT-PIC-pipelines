#!/home/xnat/anaconda3/envs/py27/bin/python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 16:30:49 2019

@author: xnat
"""

from pyxnat import Interface
import os


def upload_XNAT_files(download_info, path, xnat_folder):
    url = "/project/"
    path = os.path.expanduser("~") + "/Documents/.central.cfg"
    for j in range(0, len(download_info)):
        for experiment_dir in download_info[j]["Scan Dirs"]:
            # print("Scan dir è {}",format(experiment_dir))
            list_files = os.listdir(
                download_info[j]["Subject Folder"]
                + "/"
                + download_info[j]["Session ID"]
                + "/"
                + download_info[j]["Scan ID"]
            )
	    folder2xnat = (download_info[j]["Subject Folder"]
        			+ "/"
        			+ download_info[j]["Session ID"]
        			+ "/"
        			+ download_info[j]["Scan ID"])
	    print('folder2xnat is',folder2xnat)	
            # file_folder=(download_info[j]["Subject Folder"]+'/'+scan_dirs)
            # tmp=file_folder.split('PROCESSED/')
            # tmp=tmp[1].split('/')
            # tmp=tmp[1].split('/')
            # scan_id=tmp[0].split('-')[0]
            # print(list_files)
            # print(file_folder)
            for files in list_files:
                #file_path=download_info[j]["Subject Folder"]+'/'+scan_dirs+'/'+files
                #print(file_path)
                central = Interface(config=path)
                try:
                    serverurl = (
                        url
                        + download_info[j]["Project"]
                        + "/subject/"
                        + download_info[j]["Subject ID"]
                        + "/experiments/"
                        + download_info[j]["Session ID"]
                        + "/"
                    )
                    experiment = central.select(serverurl)
		    os.chdir(folder2xnat)
                    experiment.resource(xnat_folder).file(files).insert(files)
                finally:
                    central.disconnect()
