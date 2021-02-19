#!/bin/bash -e

if [ "$1" == "--version" ]; then
    echo 1
    exit 0
fi


input=$1
output=$2
rtstruct_file=$3
original_dcm_path=$4
mask_tumor_path=$5
T2_map_path=$6
dir=pwd

rtstruct="/RTSTRUCT"
rtstruct_folder_path="$input$rtstruct"
rtstruct_file_path=$(find $rtstruct_folder_path  -name \*.dcm)
#rtstruct_file_path="$rtstruct_folder_path$rtstruct_file"

echo -e "\n"
echo "*** run.sh ***"
echo -e "\n"
echo "RTSTRUCT file is $rtstruct_file_path"

original_dcm_path=$(find $input -type d -regextype sed -regex ".*/[0-9]")
echo "Original DICOM path is $original_dcm_path"

echo  $input > "/home/xnat/input.txt"
echo $output > "/home/xnat/output.txt"

/home/xnat/anaconda3/bin/python3 /home/xnat/pipeline/xnat-pipeline-engine/templates/scripts/mask_average/rtstruct2nii.py "$rtstruct_file_path" "$original_dcm_path" "$output"

mask_tumor_path=$(find $output -type f -name "mask_*.nii")
echo "Mask tumor path is $mask_tumor_path"

T2_map_path=$(find $input -type f -name "*T2_map.nii")
echo "T2 map path is $T2_map_path"

/home/xnat/anaconda3/bin/python3 /home/xnat/pipeline/xnat-pipeline-engine/templates/scripts/mask_average/mask_average.py "$T2_map_path" "$mask_tumor_path" "$output"

echo "Done"
