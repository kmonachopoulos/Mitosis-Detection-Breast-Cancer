%-------------------------------------------------------------------------%
%  Project       : Mitosis-Detection-Breast-Cancer                        %
%  File          : mitosis_recognition.m                                  %
%  Description   : Detect Mitotic Cells in Histological Images            %
%  Author        : Monahopoulos Konstantinos                              %
%-------------------------------------------------------------------------%

clc;clear all;close all;
warning off;
categories=5;
string_handler=[3 4 5 7 10];
string_handler_char=char('A03','A04','A05','A07','A10');
mit_cell_sub_image_cnt=0;
not_mit_cell_sub_image_cnt=0;
for cat=1:categories
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dir_mitosis = sprintf('G:\\Ç.Å.Ð ÐÁÔÑÁ ÌÅÔÁÐÔÕ×ÉÁÊÏ\\ØçöéáêÞ åðåîåñãáóßá åéêüíáò êáé video\\Project\\%s\\mitosis',string_handler_char(cat,:));
    list_of_csv = dir( fullfile(dir_mitosis,'*.csv') );
    list_of_csv_names = {list_of_csv.name}';
    
    dir_frames = sprintf('G:\\Ç.Å.Ð ÐÁÔÑÁ ÌÅÔÁÐÔÕ×ÉÁÊÏ\\ØçöéáêÞ åðåîåñãáóßá åéêüíáò êáé video\\Project\\%s\\frames\\x40',string_handler_char(cat,:));
    list_of_img = dir( fullfile(dir_frames,'*.tiff') );
    list_of_img_names = {list_of_img.name}';
    
    switch cat
        
        case 1
            %Hold 1 test image from A03 category, the rest for training
            labels=ceil((1:size(list_of_img,1))./10);
            cv= cvpartition(labels,'HoldOut',0.02);
            train_label=cv.training;
        case 2
            %Hold 1 test image from A04 category, the rest for training
            labels=ceil((1:size(list_of_img,1))./10);
            cv= cvpartition(labels,'HoldOut',0.01);
            train_label=cv.training;
        case 3
            %Hold 1 test image from A05 category, the rest for training
            labels=ceil((1:size(list_of_img,1))./10);
            cv= cvpartition(labels,'HoldOut',0.01);
            train_label=cv.training;
        case 4
            %Hold 1 test image from A07 category, the rest for training
            labels=ceil((1:size(list_of_img,1))./10);
            cv= cvpartition(labels,'HoldOut',0.02);
            train_label=cv.training;
        case 5
            %Hold 1 test image from A10 category, the rest for training
            labels=ceil((1:size(list_of_img,1))./10);
            cv= cvpartition(labels,'HoldOut',0.02);
            train_label=cv.training;
            
    end
    
    gen_counter=0; % general counter
    for i=1:2:numel(list_of_csv_names)
        gen_counter=gen_counter+1;
        string=sprintf('reading image %d of category %s',gen_counter,string_handler_char(cat,:));
        disp(string);
        error=0;
        
        not_mitotic_coor_fname=fullfile(dir_mitosis,list_of_csv_names{i+1}); %Read non mitotic cells csv filename
        mitotic_coor_fname = fullfile(dir_mitosis,list_of_csv_names{i});     %Read mitotic cells csv filename
        fname_frames = fullfile(dir_frames,list_of_img_names{gen_counter});
        
        try
            
            Histology_frames(:,:,:,gen_counter) = imread(fname_frames); %Read the histology image from the specific category
            %Histology_frames_str{gen_counter}= imread(fname_frames);
            
            mitotic_cell_coordinates{gen_counter,cat} = (csvread(mitotic_coor_fname)); %Read mitotic cells csv coordinates
            not_mitotic_cell_coordinates{gen_counter,cat} = (csvread(not_mitotic_coor_fname)); %Read non mitotic cells csv coordinates
            
        catch err
            % catch error of reading empty csv coordinates file
            error=1; % change error state conditional variable
        end
        if ( train_label(gen_counter) && error == 0)
            
            string=sprintf('Image %d selected as train image ',gen_counter);
            disp(string)
            
            % Cut a sub-image (of mitotic cell) from the specific coordinates
            
            mit_cell_sub_image= GetCell(Histology_frames(:,:,:,gen_counter)... % Call function GetCell
                ,mitotic_cell_coordinates{gen_counter,cat});
            % Cut a sub-image (of non mitotic cell) from the specific coordinates
            not_mit_cell_sub_image= GetCell(Histology_frames(:,:,:,gen_counter)... % Call function GetCell
                ,not_mitotic_cell_coordinates{gen_counter,cat});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEATURE EXTRACTION - TRAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('feature extraction..');
            
            if (~(1 && ~any(mit_cell_sub_image(:)))) % Check mit_cell_sub_image if all zeros...
                for mits=1:size(mit_cell_sub_image,3)
                    mit_cell_sub_image_cnt=mit_cell_sub_image_cnt+1;
                    % feature extraction (LBP) for EVERY mitotic sub-image (need it for train)
                    LBPsMitotic(mit_cell_sub_image_cnt,:)=lbp(mit_cell_sub_image(:,:,mits));
                    
                end
            end
            if (~(1 && ~any(not_mit_cell_sub_image(:)))) % Check not_mit_cell_sub_image if all zeros...
                for mits=1:size(not_mit_cell_sub_image,3)
                    not_mit_cell_sub_image_cnt=not_mit_cell_sub_image_cnt+1;
                    % feature extraction (LBP) for EVERY non mitotic sub-image (need it for train)
                    LBPsNotMitotic(not_mit_cell_sub_image_cnt,:)=lbp(not_mit_cell_sub_image(:,:,mits));
                end
            end
        else if (train_label(gen_counter)==0)
                
                string=sprintf('Image %d selected as test image',gen_counter);
                disp(string)
                % From every test image that is selected as train image detect and extract all the cells that contains (sub-images and coordinates)...
                switch cat
                    case 1
                        disp('Random train image from category A03')
                        
                        [all_cells_A03 all_cells_coor_A03]=GetAllCells(Histology_frames(:,:,:,gen_counter)); % Call function GetAllCells
                        mytitle= sprintf('Number of cells detected  -> %d , Train Image Category -> %s',size(all_cells_coor_A03,1)...
                            ,string_handler_char(cat,:));
                        title(mytitle)
                        disp('feature extraction..');
                        % feature extraction (LBP) for every sub-image-cell of test image
                        for cell_counter=1:size(all_cells_A03,3)
                            LBPsAllCells_A03(cell_counter,:)=  lbp(all_cells_A03(:,:,cell_counter));
                        end
                        
                    case 2
                        disp('Random train image from category A04')
                        
                        [all_cells_A04 all_cells_coor_A04]=GetAllCells(Histology_frames(:,:,:,gen_counter));% Call function GetAllCells
                        mytitle= sprintf('Number of cells detected  -> %d , Train Image Category -> %s',size(all_cells_coor_A04,1)...
                            ,string_handler_char(cat,:));
                        title(mytitle)
                        disp('feature extraction..');
                        % feature extraction (LBP) for every sub-image-cell of test image
                        for cell_counter=1:size(all_cells_A04,3)
                            LBPsAllCells_A04(cell_counter,:)=  lbp(all_cells_A04(:,:,cell_counter));
                        end
                    case 3
                        disp('Random train image from category A05')
                        
                        [all_cells_A05 all_cells_coor_A05]=GetAllCells(Histology_frames(:,:,:,gen_counter));% Call function GetAllCells
                        mytitle= sprintf('Number of cells detected  -> %d , Train Image Category -> %s',size(all_cells_coor_A05,1)...
                            ,string_handler_char(cat,:));
                        title(mytitle)
                        disp('feature extraction..');
                        % feature extraction (LBP) for every sub-image-cell of test image
                        for cell_counter=1:size(all_cells_A05,3)
                            LBPsAllCells_A05(cell_counter,:)=  lbp(all_cells_A05(:,:,cell_counter));
                        end
                    case 4
                        disp('Random train image from category A07')
                        
                        [all_cells_A07 all_cells_coor_A07]=GetAllCells(Histology_frames(:,:,:,gen_counter));% Call function GetAllCells
                        mytitle= sprintf('Number of cells detected  -> %d , Train Image Category -> %s',size(all_cells_coor_A07,1)...
                            ,string_handler_char(cat,:));
                        title(mytitle)
                        disp('feature extraction..');
                        % feature extraction (LBP) for every sub-image-cell of test image
                        for cell_counter=1:size(all_cells_A07,3)
                            LBPsAllCells_A07(cell_counter,:)=  lbp(all_cells_A07(:,:,cell_counter));
                        end
                    case 5
                        disp('Random train image from category A10')
                        
                        [all_cells_A10 all_cells_coor_A10] =GetAllCells(Histology_frames(:,:,:,gen_counter));% Call function GetAllCells
                        mytitle= sprintf('Number of cells detected  -> %d , Train Image Category -> %s',size(all_cells_coor_A10,1)...
                            ,string_handler_char(cat,:));
                        title(mytitle)
                        disp('feature extraction..');
                        % feature extraction (LBP) for every sub-image-cell of test image
                        for cell_counter=1:size(all_cells_A10,3)
                            LBPsAllCells_A10(cell_counter,:)=  lbp(all_cells_A10(:,:,cell_counter));
                        end
                end
            end
        end
        disp('........')
    end
    
    disp('Change Category - Repeat ...')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mitotic_feature_vector = mean(LBPsMitotic,1); % Make a Mitotic feature vector
Non_Mitotic_feature_vector = mean(LBPsNotMitotic,1);  % Make a non Mitotic feature vector

D=0;  % Number of detected mitosis
TP=0; % Number of true positives
FP=0; % Number of False positivies
FN=0; % False negatives
micro_meter_per_pixel=0.2455; % Aperio analysis

for cat=1:categories
    
    
    switch cat
        case 1
            
            cellsize=size(LBPsAllCells_A03,1);
            for cell_counter=1:cellsize
                
                % Euclidean distance of every cell in test image with the feature vectors
                distance_to_mitotic = pdist2(LBPsAllCells_A03(cell_counter,:),Mitotic_feature_vector,'euclidean');
                distance_to_non_mitotic = pdist2(LBPsAllCells_A03(cell_counter,:),Non_Mitotic_feature_vector,'euclidean');
                
                % Decision-making if it is a mitotic or not...
                if(distance_to_mitotic<distance_to_non_mitotic) % Predicted as mitotic
                    D=D+1;
                    if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted wrong
                        disp('prediction -> mitotic , cell -> non mitotic . empty csv')
                        FP=FP+1;
                    else  % Csv file is not empty
                        for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                            
                            % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                            % point of a ground truth mitosis according to
                            % ICPR 2014 Contest...
                            if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                    all_cells_coor_A03(cell_counter,:),'euclidean')*micro_meter_per_pixel)<8)
                                TP=TP+1;
                                disp('prediction -> mitotic , cell -> mitotic')
                            else
                                
                                % If i predicted the cell as mitotic and there are mitotic cells in the image but the coordinates differ then i predicted wrong
                                FP=FP+1;
                                disp('prediction -> mitotic , cell -> non mitotic')
                            end
                            
                        end
                    end
                    
                    % Decision-making if it is a mitotic or not...
                else if (distance_to_mitotic>distance_to_non_mitotic) % Predicted as non mitotic
                        
                        if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted right
                            disp('prediction -> non mitotic , cell -> non mitotic . empty csv')
                            
                        else  % Csv file is not empty
                            for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                                
                                % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                                % point of a ground truth mitosis according to
                                % ICPR 2014 Contest...
                                if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                        all_cells_coor_A03(cell_counter,:),'euclidean')*micro_meter_per_pixel)>8)
                                    
                                    disp('prediction -> non mitotic , cell -> non mitotic')
                                else
                                    FN=FN+1;
                                    % If i predicted the cell as non mitotic and there are mitotic cells in the image with close coordinates then i predicted wrong
                                    disp('prediction -> non mitotic , cell -> mitotic')
                                end
                            end
                        end
                    end
                end
            end
        case 2
            cellsize=size(LBPsAllCells_A04,1);
            for cell_counter=1:cellsize
                
                % Euclidean distance of every cell in test image with the feature vectors
                distance_to_mitotic = pdist2(LBPsAllCells_A04(cell_counter,:),Mitotic_feature_vector,'euclidean') ;
                distance_to_non_mitotic = pdist2(LBPsAllCells_A04(cell_counter,:),Non_Mitotic_feature_vector,'euclidean') ;
                
                % Decision-making if it is a mitotic or not...
                if(distance_to_mitotic<distance_to_non_mitotic) % Predicted as mitotic
                    D=D+1;
                    if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted wrong
                        disp('prediction -> mitotic , cell -> non mitotic . empty csv')
                        FP=FP+1;
                    else  % Csv file is not empty
                        for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                            
                            % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                            % point of a ground truth mitosis according to
                            % ICPR 2014 Contest...
                            if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                    all_cells_coor_A04(cell_counter,:),'euclidean')*micro_meter_per_pixel)<8)
                                TP=TP+1;
                                disp('prediction -> mitotic , cell -> mitotic')
                            else
                                
                                % If i predicted the cell as mitotic and there are mitotic cells in the image but the coordinates differ then i predicted wrong
                                FP=FP+1;
                                disp('prediction -> mitotic , cell -> non mitotic')
                            end
                            
                        end
                    end
                    
                    % Decision-making if it is a mitotic or not...
                else if (distance_to_mitotic>distance_to_non_mitotic) % Predicted as non mitotic
                        
                        if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted right
                            disp('prediction -> non mitotic , cell -> non mitotic . empty csv')
                            
                        else  % Csv file is not empty
                            for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                                
                                % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                                % point of a ground truth mitosis according to
                                % ICPR 2014 Contest...
                                if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                        all_cells_coor_A04(cell_counter,:),'euclidean')*micro_meter_per_pixel)>8)
                                    
                                    disp('prediction -> non mitotic , cell -> non mitotic')
                                    FP=FP+1;
                                else
                                    FN=FN+1;
                                    % If i predicted the cell as non mitotic and there are mitotic cells in the image with close coordinates then i predicted wrong
                                    disp('prediction -> non mitotic , cell -> mitotic')
                                end
                            end
                        end
                    end
                end
            end
        case 3
            cellsize=size(LBPsAllCells_A05,1);
            for cell_counter=1:cellsize
                
                % Euclidean distance of every cell in test image with the feature vectors
                distance_to_mitotic = pdist2(LBPsAllCells_A05(cell_counter,:),Mitotic_feature_vector,'euclidean') ;
                distance_to_non_mitotic = pdist2(LBPsAllCells_A05(cell_counter,:),Non_Mitotic_feature_vector,'euclidean') ;
                
                % Decision-making if it is a mitotic or not...
                if(distance_to_mitotic<distance_to_non_mitotic) % Predicted as mitotic
                    D=D+1;
                    if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted wrong
                        disp('prediction -> mitotic , cell -> non mitotic . empty csv')
                        FP=FP+1;
                    else  % Csv file is not empty
                        for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                            
                            % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                            % point of a ground truth mitosis according to
                            % ICPR 2014 Contest...
                            if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                    all_cells_coor_A05(cell_counter,:),'euclidean')*micro_meter_per_pixel)<8)
                                TP=TP+1;
                                disp('prediction -> mitotic , cell -> mitotic')
                            else
                                
                                % If i predicted the cell as mitotic and there are mitotic cells in the image but the coordinates differ then i predicted wrong
                                FP=FP+1;
                                disp('prediction -> mitotic , cell -> non mitotic')
                            end
                            
                        end
                    end
                    
                    % Decision-making if it is a mitotic or not...
                else if (distance_to_mitotic>distance_to_non_mitotic) % Predicted as non mitotic
                        
                        if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted right
                            disp('prediction -> non mitotic , cell -> non mitotic . empty csv')
                            
                        else  % Csv file is not empty
                            for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                                
                                % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                                % point of a ground truth mitosis according to
                                % ICPR 2014 Contest...
                                if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                        all_cells_coor_A05(cell_counter,:),'euclidean')*micro_meter_per_pixel)>8)
                                    
                                    disp('prediction -> non mitotic , cell -> non mitotic')
                                else
                                    FN=FN+1;
                                    % If i predicted the cell as non mitotic and there are mitotic cells in the image with close coordinates then i predicted wrong
                                    disp('prediction -> non mitotic , cell -> mitotic')
                                end
                            end
                        end
                    end
                end
            end
        case 4
            cellsize=size(LBPsAllCells_A07,1);
            for cell_counter=1:cellsize
                
                % Euclidean distance of every cell in test image with the feature vectors
                distance_to_mitotic = pdist2(LBPsAllCells_A07(cell_counter,:),Mitotic_feature_vector,'euclidean') ;
                distance_to_non_mitotic = pdist2(LBPsAllCells_A07(cell_counter,:),Non_Mitotic_feature_vector,'euclidean') ;
                
                % Decision-making if it is a mitotic or not...
                if(distance_to_mitotic<distance_to_non_mitotic) % Predicted as mitotic
                    D=D+1;
                    if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted wrong
                        disp('prediction -> mitotic , cell -> non mitotic . empty csv')
                        FP=FP+1;
                    else  % Csv file is not empty
                        for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                            
                            % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                            % point of a ground truth mitosis according to
                            % ICPR 2014 Contest...
                            if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                    all_cells_coor_A07(cell_counter,:),'euclidean')*micro_meter_per_pixel)<8)
                                TP=TP+1;
                                disp('prediction -> mitotic , cell -> mitotic')
                            else
                                
                                % If i predicted the cell as mitotic and there are mitotic cells in the image but the coordinates differ then i predicted wrong
                                
                                disp('prediction -> mitotic , cell -> non mitotic')
                            end
                            
                        end
                    end
                    
                    % Decision-making if it is a mitotic or not...
                else if (distance_to_mitotic>distance_to_non_mitotic) % Predicted as non mitotic
                        
                        if( isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted right
                            disp('prediction -> non mitotic , cell -> non mitotic . empty csv')
                            
                        else  % Csv file is not empty
                            for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                                
                                % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                                % point of a ground truth mitosis according to
                                % ICPR 2014 Contest...
                                if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                        all_cells_coor_A07(cell_counter,:),'euclidean')*micro_meter_per_pixel)>8)
                                    
                                    disp('prediction -> non mitotic , cell -> non mitotic')
                                else
                                    FN=FN+1;
                                    % If i predicted the cell as non mitotic and there are mitotic cells in the image with close coordinates then i predicted wrong
                                    disp('prediction -> non mitotic , cell -> mitotic')
                                end
                            end
                        end
                    end
                end
            end
        case 5
            cellsize=size(LBPsAllCells_A10,1);
            for cell_counter=1:cellsize
                
                % Euclidean distance of every cell in test image with the feature vectors
                distance_to_mitotic = pdist2(LBPsAllCells_A10(cell_counter,:),Mitotic_feature_vector,'euclidean') ;
                distance_to_non_mitotic = pdist2(LBPsAllCells_A10(cell_counter,:),Non_Mitotic_feature_vector,'euclidean') ;
                
                % Decision-making if it is a mitotic or not...
                if(distance_to_mitotic<distance_to_non_mitotic) % Predicted as mitotic
                    D=D+1;
                    if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted wrong
                        disp('prediction -> mitotic , cell -> non mitotic . empty csv')
                        FP=FP+1;
                    else  % Csv file is not empty
                        for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                            
                            % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                            % point of a ground truth mitosis according to
                            % ICPR 2014 Contest...
                            if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                    all_cells_coor_A10(cell_counter,:),'euclidean')*micro_meter_per_pixel)<8)
                                TP=TP+1;
                                disp('prediction -> mitotic , cell -> mitotic')
                            else
                                
                                % If i predicted the cell as mitotic and there are mitotic cells in the image but the coordinates differ then i predicted wrong
                                FP=FP+1;
                                disp('prediction -> mitotic , cell -> non mitotic')
                            end
                            
                        end
                    end
                    
                    % Decision-making if it is a mitotic or not...
                else if (distance_to_mitotic>distance_to_non_mitotic) % Predicted as non mitotic
                        
                        if(isempty(mitotic_cell_coordinates{find(train_label==0),cat})); % Empty csv file no mitotic cells in this image, i predicted right
                            disp('prediction -> non mitotic , cell -> non mitotic . empty csv')
                            
                        else  % Csv file is not empty
                            for gen_counter=1:size(mitotic_cell_coordinates{find(train_label==0),cat},1)
                                
                                % A detected mitosis would be counted as correct if its centre point is localised within a range of 8 micrometer of the centre
                                % point of a ground truth mitosis according to
                                % ICPR 2014 Contest...
                                if((pdist2(mitotic_cell_coordinates{find(train_label==0),cat}(gen_counter,1:2),...
                                        all_cells_coor_A10(cell_counter,:),'euclidean')*micro_meter_per_pixel)>8)
                                    
                                    disp('prediction -> non mitotic , cell -> non mitotic')
                                else
                                    FN=FN+1;
                                    % If i predicted the cell as non mitotic and there are mitotic cells in the image with close coordinates then i predicted wrong
                                    disp('prediction -> non mitotic , cell -> mitotic')
                                end
                            end
                        end
                    end
                end
            end
    end
end

Recall=TP/(TP+FN);
Precision=TP/(TP+FP);
Fmeasure=2*((Precision*Recall)/(Precision+Recall));
mytitle= sprintf('Recall is %.3f %%',Recall*100);
disp(mytitle)
mytitle= sprintf('Precision is %.3f %%',Precision*100);
disp(mytitle)
mytitle= sprintf('Fmeasure is %.3f %%',Fmeasure*100);
disp(mytitle)
mytitle= sprintf('F1-score is %.3f %%',((Precision+Recall)/2)*100);
disp(mytitle)