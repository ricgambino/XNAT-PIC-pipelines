% DWI analysis for DICOM files
%
% authors: Sara Zullino and Dario Longo
% 26/02/2019
% 

function process_DWI(input_img_dir,output_img_dir)

    % Extract matrix from DICOM files
    %img_dir=uigetdir('','Select DICOM folder (number) where DWI image is stored');
    [pathstr, folder_name] = fileparts(input_img_dir);
    Files=dir(fullfile(input_img_dir,'*.dcm'));
    N = natsortfiles({Files.name});
    for ii = 1:numel(Files)
       matrix(:,:,:,ii) = dicomread(fullfile(input_img_dir,N{ii})); 
    end

    image(:,:,:) = squeeze(matrix(:,:,:));

    %% Look for DICOM tags needed for processing

    size_x = size(image,1);
    size_y = size(image,2);
    nFrames = size(image,3);     % The number of files in the directory
    size_xy = [size_x size_y];
    
    factor = double((2^(32)-1)/double(max(max(max(image(:,:,:)))))); 
    
    for k=1:numel(Files)
        FileName = N{k};
        metadata = dicominfo( fullfile(fullfile(input_img_dir,FileName)));
        
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
        nExps = nFrames_dicom./nSlices;

        % b-values   
        b_val(1,k) = metadata.(dicomlookup('0018', '9087'));
        
        % Slice thickness
        slice_thickness(1,k) = metadata.(dicomlookup('0018', '0050'));

        % Slope
        slope(1,k) = metadata.(dicomlookup('0028', '1053'));
        
        % Factor
        dicom_factor = metadata.(dicomlookup('0028', '1053')); 
        
        % Image conversion from 16 to 32 bit
        image_vect(:,:,k)=double((double(image(:,:,k))./double(dicom_factor)));

    end
      
    size_z = nFrames./nExps;
    reshape_conv_img=reshape(image_vect,size_x,size_y,size_z,nExps);
    
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

    % Compute ADC
    b_values=b_val(1,1:size_z:end);
    start_point1=1.0e-3;
    start_point2=[1.0e-3 0.1];
    start_point3=[1.0e-3 0.1 1.0e-2];
    vect_analysis = zeros(1,nExps);
    ADC_image = zeros(size_x,size_y,1);

    for z=1:nSlices
        for i=1:size_y 
            for j=1:size_x
                if seg_img(i,j,z)>0
                   for t=1:nExps
                       vect_analysis(t)=reshape_conv_img(i,j,z,t)/reshape_conv_img(i,j,z,1);
                   end
                   [estimates_ADC] = fit_ADC(b_values,vect_analysis, start_point1);
                   ADC_image(i,j,z) = estimates_ADC;
                else
                   ADC_image(i,j,z)=0;
                end
                clear vect_analysis
            end
        end
    end

    cd(output_img_dir);
    
    % ADC map figure
    h_adcmap=figure('Name', 'ADC map [mm2/sec]');
    nElems=single(size_z);
    subplotRows=1;
    subplotRows=max(1, subplotRows);
    subplotCols=ceil(nElems/subplotRows);
    for iElem=1:nElems
       figure(h_adcmap);
       ax(iElem)=subplot_tight(subplotRows, subplotCols, iElem, [0.009]); % distance between images
       imshow(ADC_image(:,:,iElem),[1e-4 1e-2]);
       title(['slice ' num2str(iElem)])
       h=colorbar;
       set(h, 'Position', [.92 .301 .03 .398])  %[h pos,v pos,thickness,bar length]
       set(h,'FontSize',14, 'FontWeight', 'bold');
       pos=get(ax(iElem), 'Position');
       set(ax(iElem), 'Position', [0.89*pos(1) pos(2) pos(3) pos(4)]);
    end
    saveas(h_adcmap, [folder_name '_ADC_map.jpg']);
    
    
    % Save nifti
    resolution=[fov(1)*10/size_xy(1), fov(2)*10/size_xy(2), slice_thickness(1)];
    origin =[0 0 0]; 
    datatype =64;
    
    adc_nii= make_nii(ADC_image(:,:,:), origin, datatype, '');
    adc_nii_or=rri_orient_t1w(adc_nii);
    save_nii(adc_nii_or, [folder_name '_ADC_map.nii']);

    save([folder_name,'_DWI_IVIM_map.mat']);
    
    fileout=[folder_name '_log_file.txt']
    fid=fopen([fileout],'a+');
    %if fid<0
    %   msgbox('Problems with file name: check and try again');
    %   return
    %end

    fprintf(fid,'----------------------------------------------------------------\n');
    fprintf(fid, 'script: process DWI \n');
    fprintf(fid, 'Authors: Sara Zullino, Alessandro Paglialonga and Dario Longo\n');
    fprintf(fid,'----------------------------------------------------------------\n');
    fprintf(fid,'DWI analysis for DICOM files\n');
    fprintf(fid,'1. Extract matrix from DICOM files\n');
    fprintf(fid,'2. Look for DICOM tags needed for processing\n');
    fprintf(fid,'3. Global histogram threshold using Otsu method\n');
    fprintf(fid,'4. Compute ADC\n');
    fprintf(fid,'5. Produce output files\n');
    fclose(fid);

end
