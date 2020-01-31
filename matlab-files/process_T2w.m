%% T2 weighted map for DICOM files
%
% authors: Sara Zullino and Dario Longo
% 16/04/2019
% 
function process_T2w(input_img_dir,output_img_dir)

    % Extract matrix from DICOM files
    %input_img_dir=uigetdir('','Select DICOM folder (number) where DWI image is stored');
    [pathstr, folder_name] = fileparts(input_img_dir);

    Files=dir(fullfile(input_img_dir,'*.dcm'));
    N = natsortfiles({Files.name});
    for ii = 1:numel(Files)
       matrix(:,:,:,ii) = dicomread(fullfile(input_img_dir,N{ii})); 
    end

    image(:,:,:) = squeeze(matrix(:,:,:));

    %% Look for DICOM tags needed for processing

    size_x = size(image,1)
    size_y = size(image,2)
    nFrames = size(image,3)     % The number of files in the directory
    size_xy = [size_x size_y];
    
    factor = double((2^(32)-1)/double(max(max(max(image(:,:,:)))))); 
    
    for k=1:numel(Files)
        FileName = N{k};
        metadata = dicominfo( fullfile(fullfile(input_img_dir,FileName)));
        
        % Date
        date = metadata.(dicomlookup('0008','0022'));
        
        % Time
        time = metadata.(dicomlookup('0008','0032'));
        
        % Rows
        size_x_dicom = metadata.(dicomlookup('0028','0010')); 
        
        % Columns
        size_y_dicom = metadata.(dicomlookup('0028','0011'));     
        
        % Images in Acquisitions
        nFrames_dicom = metadata.(dicomlookup('0020','1002')); 

        if (size_x~=size_x_dicom) | (size_y~=size_y_dicom) | (nFrames~=nFrames_dicom)
            error('Matrix size does not match!');
        end
        
        % Image orientation
        img_orientation{k} = metadata.(dicomlookup('0020', '0020'));

        % Reconstruction Field of View
        fov(k,:) = metadata.(dicomlookup('0018', '9317'))';

        % Number of slices
        nSlices = metadata.(dicomlookup('0054', '0081'));
               
        % Number of experiments
        nEchoes = nFrames_dicom./nSlices;
       
        % Slice thickness
        slice_thickness(1,k) = metadata.(dicomlookup('0018', '0050'));
        
        % Techo-values   
        tEchoes_temp(1,k) = metadata.(dicomlookup('0018', '0081'));

        % Slope
        slope(1,k) = metadata.(dicomlookup('0028', '1053'));
        
        % Factor
        dicom_factor = metadata.(dicomlookup('0028', '1053')); 
        
        % Image conversion from 16 to 32 bit
        image_vect(:,:,k) = double((double(image(:,:,k))./double(dicom_factor)));

    end
    
    tEchoes = tEchoes_temp(1:nSlices:nFrames)
          
    % Global histogram threshold using Otsu's method
    seg_img = zeros(size_x,size_y,nSlices);
    for h=1:nSlices
        [counts,x] = imhist(image(:,:,h),16);
        threshold = otsuthresh(counts);
        seg_img(:,:,h) = im2bw(image(:,:,h),threshold);
        clear threshold;
        clear x;
        clear counts;
    end

    % Compute T2
    vect_analysis = zeros(1,nEchoes);
    for i=1:size_y 
        for j=1:size_x
            if seg_img(i,j,1)>0
               for t=1:nEchoes
                   vect_analysis(t) = image_vect(i,j,t);
               end
               start_point = [(image_vect(i,j,1)-image_vect(i,j,nEchoes)), 10, image_vect(i,j,nEchoes)];
               [estimates, model, exitflag] = fitcurveT2(tEchoes, vect_analysis, start_point);
               if (exitflag < 1) %iteration limit reached, not converged
                   T2w_image(i,j) = -1;
                   Mzero_image(i,j) = -1;
                   off_image(i,j) = -1;
               elseif (estimates(2) > 1000)  %stima errata
                   T2w_image(i,j) = -2;
                   Mzero_image(i,j) = -2;
                   off_image(i,j) = -2;
               else
                    T2w_image(i,j) = estimates(2);
                    Mzero_image(i,j) = estimates(1);
                    off_image(i,j) = estimates(3);
               end
            else
               T2w_image(i,j) = 0;
               Mzero_image(i,j) = 0;
               off_image(i,j) = 0;
            end
        end
    end

    cd(output_img_dir);
    %mkdir('T2_images');
    %cd('T2_images');   
    
    % T2 map figure
    h_t2w_map = figure('Name', 'T2w map [ms]');
    T2w_map = imshow(T2w_image, [0 5]);
    axis square   
    colorbar
    title('T2w map [ms]')
    saveas(h_t2w_map, [folder_name '_T2w_map.jpg']);
    
    % Save nifti
    voxel_size = [fov(1)*10/size_xy(1), fov(2)*10/size_xy(2), slice_thickness(1)];
    origin = [0 0 0]; 
    datatype = 64;
    
    T2w_nii = make_nii(T2w_image(:,:,:), voxel_size, origin, datatype, '');
    T2w_nii_or = rri_orient_no_ui(T2w_nii);
    save_nii(T2w_nii_or, [folder_name '_T2w_map.nii']);

    % R2 map figure
    for i=1:size_y 
        for j=1:size_x  
            if seg_img(i,j,1)>0 
                R2_image(i,j) = (T2w_image(i,j))^-1; 
            else
                R2_image(i,j) = 0; 
            end
        end
    end

    h_r2map = figure('Name', 'R2 map [s]');
    imshow(R2_image, [0 5]); 
    axis square;
    colorbar;
    title('R2 map [s]');
    saveas(h_r2map, [folder_name '_R2_map' '.jpg']);
    
    % Save nifti
    voxel_size = [fov(1)*10/size_xy(1), fov(2)*10/size_xy(2), slice_thickness(1)];
    origin = [0 0 0]; 
    datatype = 64;
    
    R2_nii = make_nii(R2_image(:,:,:), voxel_size, origin, datatype, '');
    R2_nii_or = rri_orient_no_ui(R2_nii);
    save_nii(R2_nii_or, [folder_name '_R2_map.nii']);

fileout = [folder_name '_log_file.txt']
fid=fopen([fileout],'a+');
if fid<0
   msgbox('Problems with file name: check and try again');
   return
end

fprintf(fid,'----------------------------------------------------------------\n');
fprintf(fid, 'CIM-XNAT Pipeline: Process T2 weighted images\n');
fprintf(fid, 'Authors: Sara Zullino, Dario Longo\n');
fprintf(fid,'----------------------------------------------------------------\n');
fprintf(fid,'T2 analysis for DICOM files\n');
fprintf(fid,'1. Extract matrix from DICOM files\n');
fprintf(fid,'2. Look for DICOM tags needed for processing\n');
fprintf(fid,'3. Global histogram threshold using Otsu method\n');
fprintf(fid,'4. Compute T2\n');
fprintf(fid,'5. Produce output files\n');
fclose(fid);

% save variables
save([folder_name '_T2w_map' '_variables.mat' ])

end






