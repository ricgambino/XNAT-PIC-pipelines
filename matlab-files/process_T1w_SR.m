%% T1 Sat Rec map for DICOM files
%
% authors: Sara Zullino and Dario Longo
% 03/07/2019
% 
function process_T1_SR(input_img_dir,output_img_dir)

    % Extract matrix from DICOM files
    % input_img_dir=uigetdir('','Select DICOM folder (number) where T1 Sat Rec image is stored');
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
        nRep = nFrames_dicom./nSlices;
       
        % Slice thickness
        slice_thickness(1,k) = metadata.(dicomlookup('0018', '0050'));
        
        % Techo-values   
        tEchoes_temp(1,k) = metadata.(dicomlookup('0018', '0081'));
        
        % Trepetition-values   
        tRep_temp(1,k) = metadata.(dicomlookup('0018', '0080'));

        % Slope
        slope(1,k) = metadata.(dicomlookup('0028', '1053'));
        
        % Factor
        dicom_factor = metadata.(dicomlookup('0028', '1053')); 
        
        % Image conversion from 16 to 32 bit
        image_vect(:,:,k)=double((double(image(:,:,k))./double(dicom_factor)));

    end
    
    tRep = tRep_temp(1:nSlices:nFrames)
          
    % Global histogram threshold using Otsu's method
    seg_img = zeros(size_x,size_y,nSlices);
    for h=1:nSlices
        % Calculate a 32-bin histogram for the image.
        [counts,x] = imhist(image(:,:,h),32); 
        threshold = otsuthresh(counts);
        seg_img(:,:,h) = im2bw(image(:,:,h),threshold);
        clear threshold;
        clear x;
        clear counts;
    end

    % Compute T1_SR
    T1_image=zeros(size_y,size_x);
    
    for i=1:size_y 
        for j=1:size_x
            for z=1:nRep
                vect_analysis(z)=image_vect(i,j,z);
            end
            %if image_vect(i,j,1)>50000
            if seg_img(i,j,1)>0
               start_point=[(image_vect(i,j,1)-image_vect(i,j,nRep)), 100, image_vect(i,j,nRep)];
               [estimates, model, exitflag] = fitcurveT1(tRep, vect_analysis, start_point);
               T1_image(i,j)=estimates(2);
            else
               T1_image(i,j)=0;
            end
            clear vect_analysis
        end
    end
	
    cd(output_img_dir);

    % T1 map figure
    h_t1map=figure('Name', 'T1 SR map [ms]');
    T1_map=imshow(T1_image, [0 4000]);
    axis square   
    colorbar
    title('T1 Sat Rec map [ms]')
    saveas(h_t1map, [folder_name '_T1_SR_map.jpg']);
    
    % Save nifti
    resolution=[fov(1)*10/size_xy(1), fov(2)*10/size_xy(2), slice_thickness(1)];
    origin =[0 0 0]; 
    datatype =64;
    
    T1_nii= make_nii(T1_image(:,:,:), origin, datatype, '');
    T1_nii_or=rri_orient_t1w(T1_nii);
    save_nii(T1_nii_or, [folder_name '_T1_SR_map.nii']);

    %figura mappa R1
    for i=1:size_y 
        for j=1:size_x  
            if seg_img(i,j,1)>0 
                R1_image(i,j)= (T1_image(i,j)/1000)^-1; 
            else
                R1_image(i,j)=0; 
            end
        end
    end

    h_r1map=figure('Name', 'R1 Sat Rec map [s]');
    imshow(R1_image, [0 200]); 
    axis square;
    colorbar;
    title('R1 Sat Rec map [s]');
    saveas(h_r1map, [folder_name '_R1_SR_map' '.jpg']);
    
    % Save nifti
    resolution=[fov(1)*10/size_xy(1), fov(2)*10/size_xy(2), slice_thickness(1)];
    origin =[0 0 0]; 
    datatype =64;
    
    R1_nii= make_nii(R1_image(:,:,:), origin, datatype, '');
    R1_nii_or=rri_orient_t1w(R1_nii);
    save_nii(R1_nii_or, [folder_name '_R1_SR_map.nii']);
    
    % Log files
    fileout=[folder_name '_log_file.txt']
    fid=fopen([fileout],'a+');
    if fid<0
       msgbox('Problems with file name: check and try again');
       return
    end

    fprintf(fid,'----------------------------------------------------------------\n');
    fprintf(fid, 'script: process T1 Sat Rec \n');
    fprintf(fid, 'Authors: Sara Zullino, Alessandro Paglialonga and Dario Longo\n');
    fprintf(fid,'----------------------------------------------------------------\n');
    fprintf(fid,'T1 Sat Rec analysis for DICOM files\n');
    fprintf(fid,'1. Extract matrix from DICOM files\n');
    fprintf(fid,'2. Look for DICOM tags needed for processing\n');
    fprintf(fid,'3. Global histogram threshold using Otsu method\n');
    fprintf(fid,'4. Compute T1\n');
    fprintf(fid,'5. Produce output files\n');
    fclose(fid);

    % save variables
    save([folder_name '_T1_SR_map' '_variables.mat' ])
    
end
