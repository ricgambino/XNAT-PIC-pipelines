#!/home/xnat/anaconda3/envs/py27/bin/python2

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 16:27:07 2019

@author: xnat
"""
import requests
import os
from zipfile import ZipFile
import shutil
from pyxnat import Interface

def REST_XNAT_Getscans_bytype(download_info, path):
    url = "http://localhost:8080/data/projects/"
    resources = (
        "/resources/DICOM/files?format=zip"
    )  # If I remove /resources/DICOM I'll get SNAPSHOTS too.
    
    for j in range(0, len(download_info)):
        
        
        for scan_type in download_info[j]["Scan Type"]:
            download_info[j]["Scan Dirs"] = []
            scans = "/scans/" + scan_type #Lo '0' serve a prenderlo come elemento e non come lista composta di 1 elemento.
            serverurl = (
                url
                + download_info[j]["Project"]
                + "/subjects/"
                + download_info[j]["Subject ID"]
                + "/experiments/"
                + download_info[j]["Session ID"]
                + scans
                + resources
            )
            
            #new_url='/project/'+download_info[j]["Project"]+'/subject/'+download_info[j]["Subject ID"]+'/experiments/'+download_info[j]["Session ID"]+scans+resources
            
            r = requests.get(serverurl, auth=("admin", "xn@t4cim2020"))
            #scans=central.get(new_url)
            
            #filename = path + "/archive" + str(j) + ".zip"
    
            if r.status_code >= 200 and r.status_code < 300:
                filename = path + "/archive" + str(j) + ".zip"
                with open(filename, "wb") as fd:
                    for chunk in r.iter_content(chunk_size=128):
                        fd.write(chunk)
                # For each Subject create a folder and then extract the corresponding DICOM files inside it
                # From now on the code simply copies data from the unzipped folder to the folder we created in order to have a simple hierarchy: 'subject/
                with ZipFile(filename, "r") as zip_ref:
                    #                 subj_folder=path+"/"+download_info[j]["Subject Label"]+"/"
                    #                 if not os.path.exists(subj_folder):
                    #                     os.mkdir(subj_folder)
                    #                     pass
                    zip_ref.extractall(path)
    
                    old_paths = []
                    new_paths = []
                    for file_paths in zip_ref.namelist():
    
                        base_path = os.path.dirname(file_paths)
    
                        # print("Base path is{}".format(base_path))
                        last_folder = os.path.basename(base_path)
                        # print("last folder is {}".format(last_folder))
                        base_path_without_last_folder = os.path.dirname(base_path)
                        # print("Base path without files {}".format(base_path_without_last_folder))
                        ## Get scan ID
                        tmp = base_path.split("scans/")
                        tmp = tmp[1].split("/")
                        tmp = tmp[0].split("/")
    
                        download_info[j]["Scan ID"] = tmp[0].split("-")[0]
    
                        old_path = os.path.join(path, base_path)
                        # new_path=path+base_path_without_last_folder+'/'+download_info[j]["Scan ID"]
                        new_path = os.path.join(
                            path, base_path_without_last_folder, download_info[j]["Scan ID"]
                        )
                        # print(new_path)
                        if old_path not in old_paths:
                            # print("sono dentro")
                            old_paths.append(old_path)
                            new_paths.append(new_path)
    
                        # Rename last folder as Scan ID
                    # print(new_paths)
                    for src_path, dst_path in zip(old_paths, new_paths):
                        os.makedirs(dst_path)
                        for file in os.listdir(src_path):
                            shutil.move(os.path.join(src_path, file), dst_path + "/")
                        # print("Source path is {}, dest path is {}".format(src_path,dst_path))
                        download_info[j]["Scan Dirs"].append(dst_path)
                        shutil.rmtree(src_path)
                        # shutil.move(src_path,dst_path)
                        # for file in os.listdir(src_path):
                        #    print(dst_path)
                        # shutil.move(os.path.join(src_path,file),dst_path)
    
                os.remove(filename)
            else:
                exit("Error" + str(r.status_code))
