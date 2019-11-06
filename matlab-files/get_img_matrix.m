function get_img_matrix(input_img_dir,output_img_dir)
    % Extract matrix from DICOM files
    %input_img_dir=uigetdir('','Select DICOM folder (number) where DWI image is stored');
    [path, folder_name] = fileparts(input_img_dir);
    	
    Files=dir(fullfile(input_img_dir,'*.dcm'));
    N = natsortfiles({Files.name});
    for ii = 1:numel(Files)
       Img(:,:,:,ii) = dicomread(fullfile(input_img_dir,N{ii})); 
    end

    matrix(:,:,:) = squeeze(Img(:,:,:));
    cd(path);
    save([output_img_dir,'_matrix.mat']);
end
