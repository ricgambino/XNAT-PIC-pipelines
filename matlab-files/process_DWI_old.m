%% DWI analysis for DICOM files
%
% authors: Sara Zullino and Dario Longo
% 28/11/2018
% 
function process_DWI(input_img_dir,output_img_dir)
    % Extract matrix from DICOM files
    %input_img_dir=uigetdir('','Select DICOM folder (number) where DWI image is stored');
    [~, folder_name] = fileparts(input_img_dir);

    Files=dir(fullfile(input_img_dir,'*.dcm'));
    N = natsortfiles({Files.name});
    for ii = 1:numel(Files)
       Img(:,:,:,ii) = dicomread(fullfile(input_img_dir,N{ii})); 
    end

    matrix(:,:,:) = squeeze(Img(:,:,:));

    %% Find matrix dimension and compare with DICOM tags

    size_x = size(matrix,1);
    size_y = size(matrix,2);
    size_slices = size(matrix,3);
    size_xy = [size_x size_y];

    %% Image segmentation

    for h=1:size_slices
        mean_mat(:,h) = mean(mean(matrix(:,:,h)));
        for i=1:size_x
            for j=1:size_y
                if matrix(i,j,h) >=  mean_mat(:,h)
                    image_seg(i,j,h) = 1;
                else 
                    image_seg(i,j,h) = 0;
                end
            end
        end
    end

    for k=1:numel(Files)
        FileName = N{k};
        metadata = dicominfo( fullfile(fullfile(input_img_dir,FileName)));

        size_x_dicom = metadata.(dicomlookup('0028','0010'));      % Rows
        size_y_dicom = metadata.(dicomlookup('0028','0011'));      % Columns
        size_slices_dicom = metadata.(dicomlookup('0020','1002')); % Images in Acquisitions

        if (size_x~=size_x_dicom) && (size_y~=size_y_dicom) && (size_slices~=size_slices_dicom)
            error('Matrix size does not match!');
        end

        img_orientation = metadata.(dicomlookup('0020', '0020'));

        % Get Reconstruction Field of View

        fov = metadata.(dicomlookup('0018', '9317'))';

        % Get Number of experiments
        nslices = metadata.(dicomlookup('0054', '0081'));
        nimages = metadata.(dicomlookup('0020', '1002'));

        nexp = nimages./nslices;

        % Get b-values   
        b_val(1,k) = metadata.(dicomlookup('0018', '9087'));

        slice_thickness = metadata.(dicomlookup('0018', '0050'));

        slope = metadata.(dicomlookup('0028', '1053'));

        factor(1,k) = double((2^(32)-1)/max(max(matrix(:,:,k))));

        image_vect(:,:,k)=(double(matrix(:,:,k)).*double(factor(1,k)));


    end


    size_z=size_slices./nexp;
    reshape_image_vect=reshape(image_vect,size_x,size_y,size_z,nexp);

    %calcola ADC

    b_values=b_val(1,1:size_z:end);
    %id_t2map = waitbar(0,'Calculating ADC map...');
    for z=3:3
        for i=1:size_y 
            for j=1:size_x
     %           waitbar(z/z);

                if image_seg(i,j,z)>0
                   for t=1:nexp
                       vect_analysis(t)=reshape_image_vect(i,j,z,t)/reshape_image_vect(i,j,z,1);
                   end

                   start_point1=1.0e-3;
                   start_point2=[1.0e-3 0.1];
                   start_point3=[1.0e-3 0.1 1.0e-2];
                   [estimates_ADC] = fit_ADC(b_values,vect_analysis, start_point1);

                   ADC_image(i,j,z) = estimates_ADC;

                else
                   ADC_image(i,j,z)=0;

                end
                clear vect_analysis
            end
        end
    end
    %close(id_t2map);

    cd(output_img_dir);

    %figura mappa ADC

    for ii=3:3
        h_adcmap(ii)=figure('visible', 'off');
        imshow(squeeze(ADC_image(:,:,ii)),[1e-4 1e-2]);
        title(['ADC map [mm2/sec] - slice ' num2str(ii)]);
        h=colorbar;
        set(h,'FontSize',14, 'FontWeight', 'bold');
        saveas(h_adcmap(ii), [folder_name '_ADC_map - slice ', num2str(ii) '.jpg']);
        saveas(h_adcmap(ii), [folder_name '_ADC_map - slice ', num2str(ii) '.fig']);
    end
    save([folder_name,'DWI_IVIM_map.mat']);

end
