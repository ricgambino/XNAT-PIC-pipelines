#!/home/xnat/anaconda3/envs/py27/bin/python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 16:28:27 2019

@author: Alessandro Paglialonga
"""
from pyxnat import Interface
import os
import re


def download_project_scans_of_type(build_dir, output_dir, scan_type):

    # PARSE Project ID.
    splitted = build_dir.split("build/")
    splitted = splitted[1].split("/")
    prj_ID = splitted[0]
    # Select the data with specific constraints
    path = os.path.expanduser("~") + "/Documents/.central.cfg"
    central = Interface(config=path)
    try:
        constraints = [("xnat:mrSessionData/project", "=", prj_ID)]
        data = central.select(
            "xnat:mrSessionData",
            [
                "xnat:mrSessionData/MR_SCAN_COUNT_AGG",
                "xnat:mrSessionData/SESSION_ID",
                "xnat:mrSessionData/SUBJECT_ID",
                "xnat:mrSessionData/SUBJECT_LABEL",
            ],
        ).where(constraints)
        
        # Load the useful data in a python dictioanry
        
        dictlist = [dict()]
        j = 0
        base_dir = output_dir

        ### Download info rinominalo come scan_payload o something similar.
        reg_expr=[]
        array_scan_types=[]
        for elem in scan_type:
            pattern="(?i)"+elem+"[^\(\)]*" #match each substring containing: scan_type followed by any character excluding '(' and ')' zero or more times.
            reg_expr.append(re.compile(pattern))
            
       
        #central.array.experiments(project_id=prj_ID,experiment_type='xnat:mrSessionData',constraints={'xnat:mrSessionData/MR_SCAN_COUNT_AGG':})    
    
        # range esclude il secondo estremo
        for i in range(0, len(data)):
            array_scan_types=[]
            for regular_expression in reg_expr: 
                scan_type_retrieved=regular_expression.findall(data[i]["mr_scan_count_agg"])
                if len(scan_type_retrieved)!=0:
                    array_scan_types+=scan_type_retrieved   
            dictlist[j]["Scan Type"] = array_scan_types
            dictlist[j]["Project"] = prj_ID
            dictlist[j]["Subject ID"] = data[i]["subject_id"]
            dictlist[j]["Session ID"] = data[i]["session_id"]
            dictlist[j]["Subject Label"] = data[i]["subject_label"]
            dictlist[j]["Subject Folder"] = (
                base_dir + "/" + dictlist[j]["Subject Label"]
            )
            j = j + 1
            dictlist.insert(j, dict())
        dictlist.pop(j)  # pop the last empty element
        return dictlist
    except Exception as err:
        print(err)
    finally:
        central.disconnect()
