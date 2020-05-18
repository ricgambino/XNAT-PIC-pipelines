#!/bin/bash -e


if [ "$1" == "--version" ]; then
    echo 1
    exit 0
fi

input=$1
output=$2
dir=pwd


echo $input > "/home/xnat/input.txt"
echo $output > "/home/xnat/output.txt"
#/usr/local/MATLAB/R2015B/bin/matlab

/usr/local/bin/matlab -nosplash -nodisplay << END || die "Matlab processing failed" 
addpath(genpath('/home/xnat/pipeline/xnat-pipeline-engine/templates/scripts/process_T2w')); 
addpath(genpath('/home/xnat/Documents/Executables/Matlab-files'));
process_T2w('$input','$output') 
exit; 
END

echo "Done"

