%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%This code removes spikes from an image. It generates the image without
%%match filtering. It works in the following way.
%%1) Generates a fundamental or ultraharmonic image depending on the user
%%selection
%%2) Crops out the unwanted middle region
%%3) Replaces spikes with the average of the background region
%%4) Reproduces an image with spikes removed and with the unwanted middle
%%region replaced with a matrix of zeros. Will later be modified to replace
%%the unwanted region with a matrix of values relative to the largest value
%%in the image matrix.
%%
%%N/B: Ensure you have the file DP310_Ch1_1.data16 in your workspace. If
%%you have another file, please change line 21 to suit your preference.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all;
close all;

%load data
fid = fopen('DP310 Ch1_1.data16','rb');
Row = 11395;
Col = 363;
rf = fread(fid, [Row, Col],'short');

fclose(fid);

Fs = 357e6; %sampling frequency
F_low = 42e6;
F_high = 48e6;

N_ord = 12; %order of filter
N_pts = 2800; %Number of useful time points
%plot(rf(:,1)); 

start_pos_phase = 6;
start_neg_phase = 7146;

%Pulse inversion
%Initialize variables rf_pos_phase, rf_neg_phase and PI_rf for speed
rf_pos_phase = zeros(N_pts+1, Col);
rf_neg_phase = zeros(N_pts+1, Col);
PI_rf = zeros(N_pts+1, Col);

for i = 1:Col
    rf_pos_phase(:,i) = rf(start_pos_phase:start_pos_phase+N_pts, i);
    rf_neg_phase(:,i) = rf(start_neg_phase:start_neg_phase+N_pts,i);
    PI_rf(:,i) = rf_pos_phase(:,i) + 1*rf_neg_phase(:,i);
end

[b, a] = butter(N_ord, [F_low/(Fs/2)], 'high');
[b1, a1] = butter(N_ord, [F_high/(Fs/2)], 'low');

%Initialize the variables rf_filt, rf_filt_only and rf_filt_only1 for speed
rf_filt = zeros(N_pts+1, Col);
rf_filt_only1 = zeros(N_pts+1, Col);
rf_filt_only = zeros(N_pts+1, Col);

for i = 1:Col
    rf_filt(:,i)=filter(b,a,PI_rf(:,i));
    rf_filt_only1(:,i)=filter(b,a,rf_neg_phase(:,i));
    rf_filt(:,i)=filter(b1,a1,rf_filt(:,i));
    rf_filt_only(:,i)=filter(b1,a1,rf_filt_only1(:,i));
end

%create a variable rf_pad, to pad the center with zeros. The new
%rf_filt is the filtered data, along with an array of zeros.
rf_pad = zeros(start_pos_phase,Col);
rf_filt=[rf_pad' rf_filt'];

%pad rf_pos_phase data with zeros, to create a new array rf_pos_phase1.
rf_pos_phase1=[rf_pad' rf_pos_phase'];
%rf_pos_phase=[rf_pad' rf_pos_phase'];
    
%generate the fundamental and the ultraharmonic images.
maximum_B_mode=max(rf_filt(:));
maximum_B_modeF=max(rf_neg_phase(:));
%Polar2cart1(rf_neg_phase,Fs);

Polar2cart1(rf_filt',Fs,maximum_B_mode);
caxis([-40 -5]);
title('Ultraharmonic Image')
figure;
Polar2cart1(rf_neg_phase,Fs,maximum_B_mode);
caxis([-60 4]);
title('Fundamental Image')
%normalize data for processing
rf_filt_new = rf_filt'; %only few rows padded with zeros.

%normalization
rf_filtered = abs(hilbert(rf_filt_new));
rf_filtered = rf_filtered/max(rf_filtered(:));

%Prepare fundamental matrix for spike removal
rf_filt_new_fund = rf_neg_phase;

%normalization
rf_filtered_fund = abs(hilbert(rf_filt_new_fund));
rf_filtered_fund = rf_filtered_fund/max(rf_filtered_fund(:));

%Get user preference
disp(' ');
disp('Please select your preference based on the following options');
disp('1) Fundamental Image spike removal');
disp('2) Ultraharmonic image spike removal');
disp('3) Fundamental and Ultraharmonic image spike removal');
disp(' ');
user_preference = input('Enter 1 for Fundamental, 2 for Ultraharmonic, 3 for both Fundamental and Ultraharmonic: ');

switch user_preference
    case 1 %Fundamental image processing
       
        number_of_roi = input('Enter the number of regions of interest: ');
        if(number_of_roi == false)
            no_ROI_fundamental(rf_filt_new, rf_neg_phase, start_pos_phase, Fs, Col);
            return
        end
        
        T_fund = 20*log10(rf_filtered_fund);
        for counter = 1:number_of_roi
            %extract unwanted region
            %figure('units','normalized','outerposition',[0 0 1 1])
            %I = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
            %colormap(gray);
            
            switch counter
                case 1
                    handle1 = figure('units','normalized','outerposition',[0 0 1 1]);
                    I = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
                    colormap(gray);
                    title('Now select the unwanted region');
                    pause
                    BW4 = roipoly;
                    [J4, K4] = find(BW4 > 0); %find the coordinates of the unwanted region
                    interest4 = [J4 K4];
                    %Remove the unwanted region
                    unwanted_row = max(interest4(:,1));
                    rf_filtered = rf_filtered(unwanted_row+1:end, :);
                    rf_filtered_fund = rf_filtered_fund(unwanted_row+1:end, :);
                %otherwise
                    %continue
            end
            
            %Fundamental spike removal process
            close(handle1);
            handle2 = figure('units','normalized','outerposition',[0 0 1 1])
            I_new_fund = imagesc(20*log10(rf_filtered_fund)); caxis([-60 4]);
            colormap(gray);
            title('Use cursors to select a spike. For the fundamental, select any spike');
            [spike_col_fund, spike_row_fund] = ginput(1);
            
            title('Now select the ROI');
            pause
            BW_fund = roipoly;
            [J_fund, K_fund] = find(BW_fund > 0);
            interest_fund = [J_fund, K_fund]; %region of interest
            
            title('Now select the background region. Be sure not to select a spiked region');
            pause
            BW3_fund = roipoly;
            [J3_fund, K3_fund] = find(BW3_fund > 0); %Find the coordinates of the background region
            interest3_fund = [J3_fund, K3_fund];
            
            %get values of region of interest
            for i = 1:size(interest3_fund, 1)
                for j = 1:size(interest3_fund, 2)
                    interest3_fund(i,3) = rf_filtered_fund(interest3_fund(i,1), interest3_fund(i,2));
                end
            end
            thresh_fund =  max(interest3_fund(:,3)); %threshold value of the background region
            
            %evaluate the spike relative to the background region
            handle3 = figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            plot(rf_filtered_fund(:,round(spike_col_fund)));
            string = num2str(round(spike_row_fund));
            string = strcat('Use the cursors to select the threshold value relative to the spike value. Your spike value is on x-value = ', string);
            title(string);
            
            [x_val_fund, y_val_fund] = ginput(1);
            thresh_roi_fund = y_val_fund;
            
            [row_fund, col_fund] = find(rf_filtered_fund > thresh_fund);
            
            %roi threshold selection
            [row_new_fund, col_new_fund] = find(rf_filtered_fund > thresh_roi_fund);
            spikes_roi_fund = [row_new_fund, col_new_fund];
            
            spikes_fund = [row_fund, col_fund];
            x_final_fund = interest_fund(:,1);
            y_final_fund = interest_fund(:,2);
            final_fund = [x_final_fund, y_final_fund];
            
            %Use set diff to find the difference. Isolate the spikes
            spikes_2_fund = setdiff(spikes_fund, final_fund, 'rows');
            
            rf_filt_new2_fund = rf_filt_new_fund(unwanted_row+1:end, :);
            
            %Remove fundamental image spikes
            for j = 1:size(spikes_2_fund, 1)
                if(spikes_2_fund(j,2) - 1 == 0)
                    column_before_fund = Col;
                else
                    column_before_fund = spikes_2_fund(j,2) - 1;
                end
                
                if(spikes_2_fund(j,2)+1 > Col)
                    column_after_fund = 1;
                else
                    column_after_fund = spikes_2_fund(j,2) + 1;
                end
                
                vect_sum_fund = rf_filt_new2_fund(spikes_2_fund(j,1), column_before_fund) + rf_filt_new2_fund(spikes_2_fund(j,1), column_after_fund);
                vect_avg_fund = vect_sum_fund/2;
                rf_filt_new2_fund(spikes_2_fund(j,1), spikes_2_fund(j,2)) = vect_avg_fund; %assign averaged value
            end
            
            rf_filtered_new_fund = abs(hilbert(rf_filt_new2_fund));
            rf_filtered_new_fund = rf_filtered_new_fund/max(rf_filtered_new_fund(:));
            close(handle2);
            close(handle3);
            handle4 = figure('units','normalized','outerposition', [0 0 1 1])
            interest_new_fund = [J_fund, K_fund];
            T_new_fund = 20*log10(rf_filtered_new_fund); 
            
            I_new_fund = imagesc(20*log10(rf_filtered_new_fund)); caxis([-60 4]);
            colormap(gray);
            title('Select background for CTR computation. Be sure to select a region larger than or equal to the size of the ROI');
            pause
            BW_new_fund = roipoly;
            [J_new_fund, K_new_fund] = find(BW_new_fund > 0);
            interest_new_2_fund = [J_new_fund, K_new_fund]; %fundamental image CTR computation
            
            for i = 1:length(J_fund)
                %fundamental
                interest_new_fund(i,3) = T_new_fund(interest_new_fund(i,1), interest_new_fund(i,2));
                interest_new_fund(i,4) = T_new_fund(interest_new_2_fund(i,1), interest_new_2_fund(i,2));
            end
            plaque_fund = mean(interest_new_fund(:,3));
            background_fund =  mean(interest_new_fund(:,4));
            
            CTR_fund = plaque_fund - background_fund;
            a_fund = [plaque_fund background_fund CTR_fund];
            
            %recreate rf_filt_new_fund
            rf_filt_new3_fund = zeros(unwanted_row, Col);
            rf_filt_new_fund = [rf_filt_new3_fund; rf_filt_new2_fund];
            
            rf_filt_new_fund = rf_filt_new_fund(1:end, :);
            
            %generate final fundamental image
            rf_pad_fund = zeros(start_pos_phase, Col);
            rf_filt_new_pad_fund = [rf_pad_fund; rf_filt_new_fund];
            rf_filt_new_pad_fund = rf_filt_new_pad_fund./max(rf_filt_new_pad_fund(:));
            
            size(rf_filt_new_pad_fund)
            close(handle4);
            figure; Polar2cart1(rf_filt_new_pad_fund, Fs, 1);
            caxis([-60 4])
            title('Fundamental image');
            disp('**********************************************');
            disp('FUNDAMENTAL IMAGE DATA');
            disp('**********************************************');
            disp(strcat('CTR:', num2str(CTR_fund)));
            disp(strcat('Background:',num2str(background_fund)));
            disp(strcat('ROI;',num2str(plaque_fund)));
            disp('**********************************************');
        end
    case 2 %Ultraharmonic image spike removal
        number_of_roi = input('Enter the number of regions of interest: ');
        if(number_of_roi == false)
            %[spikes_matrix, new_rf_filtered, old_rf_filtered] = no_ROI_ultraharmonic(rf_filt_new, Col, start_pos_phase,Fs);
            [spikes_matrix, rf_filt_new3,rf_filt_new2]= no_ROI_ultraharmonic(rf_filt_new, start_pos_phase, Fs, Col);
            return
        end
        T = 20*log10(rf_filtered);
        for counter = 1:number_of_roi
            %extract unwanted region
            %figure('units','normalized','outerposition',[0 0 1 1])
            %I = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
            %colormap(gray);
            
            switch counter
                case 1
                    handle1 = figure('units','normalized','outerposition',[0 0 1 1])
                    I = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
                    colormap(gray);
                    title('Now select the unwanted region');
                    pause
                    BW4 = roipoly;
                    [J4, K4] = find(BW4 > 0); %find the coordinates of the unwanted region
                    interest4 = [J4 K4];
                    %Remove unwanted region
                    unwanted_row = max(interest4(:,1));
                    rf_filtered = rf_filtered(unwanted_row+1:end, :);
                    rf_filtered_fund = rf_filtered_fund(unwanted_row+1:end, :);
                %otherwise
                    %continue
            end
            
            %spike selection
            close (handle1); %not needed 
            handle2 = figure('units','normalized','outerposition',[0 0 1 1])
            I_new = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
            colormap(gray);
            title('Use cursors to select spike that passes through ROI');
            [spike_col, spike_row] = ginput(1);
            
            title('Now select the ROI');
            pause
            BW = roipoly;
            [J, K] = find(BW > 0); %find the coordinates of the region of interest
            interest = [J K]; %Region of interest
            
            title('Now select the background region. Be sure not to select a spiked region');
            pause
            BW3 = roipoly;
            [J3, K3] = find(BW3 > 0); %find the coordinates of the background region
            interest3 = [J3 K3]; %background region
            
            %get values for the region of interest
            for i = 1:size(interest3,1)
                for j = 1:size(interest3,2)
                    interest3(i,3) = rf_filtered(interest3(i,1), interest3(i,2));
                end
            end
            thresh = max(interest3(:,3)); %threshold value fo the background region
            
            %evaluate the spike relative to the bacground region
            handle3 = figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            plot(rf_filtered(:,round(spike_col)));
            string = num2str(round(spike_row));
            string = strcat('Use the cursors to select the threshold value relative to the spike value. Your spike is on x-value = ',string);
            title(string);
            [x_val, y_val] = ginput(1);
            thresh_roi = y_val;
            
            [row, col] = find(rf_filtered > thresh);
            
            %roi threshold selection
            [row_new, col_new] = find(rf_filtered > thresh_roi);
            spikes_roi = [row_new, col_new];
            
            spikes = [row col];
            x_final = interest(:,1);
            y_final = interest(:,2);
            final = [x_final, y_final];
            
            %Use setdiff to find the difference. Isolates the spikes
            spikes_2 = setdiff(spikes, final, 'rows');
            
            rf_filt_new2 = rf_filt_new(unwanted_row+1:end,:);
            rf_filt_new3 = rf_filt_new(unwanted_row+1:end,:);
            
            %Remove spikes
            for i = 1:size(spikes_2, 1)
                if(spikes_2(i,2) - 1 == 0)
                    column_before = Col;
                else
                    column_before = spikes_2(i,2) - 1;
                end
                if(spikes_2(i,2)+1 > Col)
                    column_after = 1;
                else
                    column_after = spikes_2(i,2) + 1;
                end
                
                vect_sum = rf_filt_new2(spikes_2(i,1), column_before) + rf_filt_new2(spikes_2(i,1), column_after);
                
                vect_avg = vect_sum/2;
                
                rf_filt_new2(spikes_2(i,1), spikes_2(i,2)) = vect_avg;
            end
            rf_filtered_new = abs(hilbert(rf_filt_new2));
            rf_filtered_new = rf_filtered_new/max(rf_filtered_new(:));
            close(handle2);
            close(handle3);
            handle4 = figure('units','normalized','outerposition',[0 0 1 1])
            interest_new = [J K];
            
            T_new = 20*log10(rf_filtered_new);
            
            I_new = imagesc(20*log10(rf_filtered_new)); caxis([-60 -10]);
            colormap(gray);
            title('Select background for CTR computation. Be sure to select a region larger than or equal to the size of the ROI');
            pause
            BW_new = roipoly;
            [J_new, K_new] = find(BW_new > 0);
            interest_new_2 = [J_new, K_new];
            
            for i = 1:length(J)
                interest_new(i,3) = T_new(interest_new(i,1), interest_new(i,2));
                interest_new(i,4) = T_new(interest_new_2(i,1), interest_new_2(i,2));
            end
            
            plaque = mean(interest_new(:,3));
            background = mean(interest_new(:,4)); 
            
            CTR = plaque - background;
            
            a = [plaque background CTR];
            
        end
        
        %Recreate rf_filt_new
        rf_filt_new3 = zeros(unwanted_row, Col);
        rf_filt_new = [rf_filt_new3; rf_filt_new2];
        
        rf_filt_new = rf_filt_new(1:end, :);
        
        %generate final image
        rf_pad = zeros(start_pos_phase, Col);
        rf_filt_new_pad = [rf_pad' rf_filt_new'];
        rf_filt_new_pad = rf_filt_new_pad./max(rf_filt_new_pad(:));
        close(handle4);
        figure; Polar2cart1(rf_filt_new_pad', Fs, 1);
        caxis([-40 0]); title('Ultraharmonic image');
        disp('**********************************************');
        disp('ULTRAHARMONIC IMAGE DATA');
        disp('**********************************************');
        disp(strcat('CTR:',num2str(CTR)));
        disp(strcat('Background:', num2str(background)));
        disp(strcat('ROI:',num2str(plaque)));
        disp('**********************************************');
    case 3 %Both ultraharmonic and fundamental
        number_of_roi = input('Enter the number of regions of interest: ');
        T = 20*log10(rf_filtered);
        T_fund = 20*log10(rf_filtered_fund); 
        if(number_of_roi == false)
            no_ROI_ultraharmonic(rf_filt_new, start_pos_phase, Fs,Col);
            no_ROI_fundamental(rf_filt_new, rf_neg_phase, start_pos_phase, Fs, Col);
            return
        end
        for counter = 1:number_of_roi
            %extract unwanted region
            %figure('units','normalized','outerposition', [0 0 1 1])
            %I = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
            %colormap(gray);
    
            switch counter
                case 1
                    handle1 = figure('units','normalized','outerposition',[0 0 1 1])
                    I = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
                    colormap(gray);
                    title('Now select the unwanted region');
                    pause
                    BW4 = roipoly;
                    [J4, K4] = find(BW4 > 0); %find the coordinates of the unwanted region
                    interest4 = [J4 K4];
                    %Remove unwanted region
                    unwanted_row = max(interest4(:,1)); 
                    rf_filtered = rf_filtered(unwanted_row+1:end, :);
                    rf_filtered_fund = rf_filtered_fund(unwanted_row+1:end, :);
                %otherwise
                    %continue
            end
    
             %fundamental spike selection. Comment this out if you don't want to 
            %process the fundamental image
    
            close(handle1);
            handle2 = figure('units','normalized','outerposition',[0 0 1 1])
            I_new_fund = imagesc(20*log10(rf_filtered_fund)); caxis([-60 4]);
            colormap(gray);
            title('Use cursors to select a spike. For the fundamental, select any spike');
            [spike_col_fund, spike_row_fund] = ginput(1);
    
            title('Now select the ROI');
            pause
            BW_fund = roipoly;
            [J_fund, K_fund] = find(BW_fund > 0);
            interest_fund = [J_fund, K_fund]; %region of interest
    
            title('Now select the background region. Be sure not to select a spiked region');
            pause
            BW3_fund = roipoly;
            [J3_fund, K3_fund] = find(BW3_fund > 0); %Find the coordinates of the background region
            interest3_fund = [J3_fund K3_fund];
    
            %get values of region of interest
            for i = 1:size(interest3_fund,1)
                for j = 1:size(interest3_fund, 2)
                    interest3_fund(i,3) = rf_filtered_fund(interest3_fund(i,1), interest3_fund(i,2));
                end
            end
            thresh_fund = max(interest3_fund(:,3)); %threshold value of the background region
    
            %evaluate the spike relative to the background region
            handle3 = figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            plot(rf_filtered_fund(:,round(spike_col_fund)));
            string = num2str(round(spike_row_fund));
            string = strcat('Use the cursors to select the threshold value relative to the spike value. Your spike value on x-value ' ,string);
            title(string);
    
            [x_val_fund, y_val_fund] = ginput(1);
            thresh_roi_fund = y_val_fund;
    
            [row_fund, col_fund] = find(rf_filtered_fund > thresh_fund);
   
            %roi threshold selection
            [row_new_fund, col_new_fund] = find(rf_filtered_fund > thresh_roi_fund);
            spikes_roi_fund = [row_new_fund, col_new_fund];
    
            spikes_fund = [row_fund col_fund];
            x_final_fund = interest_fund(:,1);
            y_final_fund = interest_fund(:,2);
            final_fund = [x_final_fund, y_final_fund];
    
            %Use set diff to find the difference. Isolates the spikes
            spikes_2_fund = setdiff(spikes_fund, final_fund, 'rows');
    
            rf_filt_new2_fund = rf_filt_new_fund(unwanted_row+1:end, :);
    
            %Remove fundamental image spikes
            for j = 1:size(spikes_2_fund,1)
                if(spikes_2_fund(j,2) - 1 == 0)
                    column_before_fund = Col;
                else
                    column_before_fund = spikes_2_fund(j,2) - 1;
                end
        
                if(spikes_2_fund(j,2)+1 > Col)
                    column_after_fund = 1;
                else
                    column_after_fund = spikes_2_fund(j,2) + 1;
                end
        
                vect_sum_fund = rf_filt_new2_fund(spikes_2_fund(j,1), column_before_fund) + rf_filt_new2_fund(spikes_2_fund(j,1), column_after_fund);
                vect_avg_fund = vect_sum_fund/2;
                rf_filt_new2_fund(spikes_2_fund(j,1), spikes_2_fund(j,2)) = vect_avg_fund;
        
            end
    
            %%same for fundamental
            rf_filtered_new_fund = abs(hilbert(rf_filt_new2_fund));
            rf_filtered_new_fund = rf_filtered_new_fund/max(rf_filtered_new_fund(:));
    
            close(handle2);
            close(handle3);
            handle4 = figure('units','normalized','outerposition',[0 0 1 1])
            interest_new_fund = [J_fund K_fund];
            T_new_fund = 20*log10(rf_filtered_new_fund);
        
            I_new_fund = imagesc(20*log10(rf_filtered_new_fund)); caxis([-60 4]);
            colormap(gray);
            title('Select background for CTR computation. Be sure to select a region larger than or equal to the size of the ROI');
            pause
            BW_new_fund = roipoly;
            [J_new_fund, K_new_fund] = find(BW_new_fund > 0);
            interest_new_2_fund = [J_new_fund K_new_fund]; %fundamental image CTR computation??
    
            for i = 1:length(J_fund)
                %fundamental
                interest_new_fund(i,3) = T_new_fund(interest_new_fund(i,1), interest_new_fund(i,2));
                interest_new_fund(i,4) = T_new_fund(interest_new_2_fund(i,1), interest_new_2_fund(i,2));
            end
    
            plaque_fund = mean(interest_new_fund(:,3));
            background_fund = mean(interest_new_fund(:,4));
    
            CTR_fund = plaque_fund - background_fund;
            a_fund = [plaque_fund background_fund CTR_fund];
    
            %recreate rf_filt_new_fund
            rf_filt_new3_fund = zeros(unwanted_row, Col);
            rf_filt_new_fund = [rf_filt_new3_fund; rf_filt_new2_fund];

            rf_filt_new_fund = rf_filt_new_fund(1:end, :);
    
            %generate final fundamental image
            rf_pad_fund = zeros(start_pos_phase, Col);
            rf_filt_new_pad_fund = [rf_pad_fund; rf_filt_new_fund];
            rf_filt_new_pad_fund = rf_filt_new_pad_fund./max(rf_filt_new_pad_fund(:));

            size(rf_filt_new_pad_fund)
            close(handle4);
            figure; Polar2cart1(rf_filt_new_pad_fund, Fs, 1);
            caxis([-60 4]); title('Fundamental image');
      
            %end of fundamental image code
    
            %return
    
            %spike selection
            handle5 = figure('units','normalized','outerposition',[0 0 1 1])
            I_new = imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
            colormap(gray);
            title('Use cursors to select spike that passes through ROI');
            [spike_col, spike_row] = ginput(1);
    
            title('Now select the ROI');
            pause
            BW = roipoly;
            [J, K] = find(BW > 0); %find the coordinates of the region of interest
            interest = [J K]; %region of interest 
            %interest_fund = [J K];
    
            title('Now select the background region. Be sure not to select a spiked region');
            pause
            BW3 = roipoly;
            [J3, K3] = find(BW3 > 0); %find the coordinates of the background region
            interest3 = [J3 K3]; %background region
            %interest3_fund = [J3 K3];
    
            %get values of region of interest
            for i = 1:size(interest3, 1)
                for j = 1:size(interest3,2)
                    interest3(i,3) = rf_filtered(interest3(i,1), interest3(i,2));
                %interest3_fund(i,3) = rf_filtered_fund(interest3_fund(i,1), interest3_fund(i,2));
            
                end
            end
            thresh = max(interest3(:,3)); %threshold value of the background region
            %thresh_fund = max(interest3_fund(:,3)); %threshold value of background region of fundamental matrix
    
            %evaluete the spike relative to the background region
            handle6 = figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            plot(rf_filtered(:,round(spike_col)));
            string2 = num2str(round(spike_row));
            string2 = strcat('Use the cursors to selct the threshold value relative to the spike value. Your spike x-value = ',string2);
            title(string2);
            [x_val, y_val] = ginput(1);
            thresh_roi = y_val;
    
            [row, col] = find(rf_filtered > thresh);
            %[row_fund col_fund] = find(rf_filtered_fund > thresh_fund);
    
            %roi threshold selection
            [row_new, col_new] = find(rf_filtered > thresh_roi);
            spikes_roi = [row_new, col_new];
    
            spikes = [row col];
            x_final = interest(:,1);
            y_final = interest(:,2);
            final = [x_final, y_final];
   
    
            %fundamental spike- evaluate the spike relative to the background
            %region
            %figure('units','normalized','outerposition',[0 0 1 1])
            %hold on
            %plot(rf_filtered_fund(:, round(spike_col)));
            %title('Use the cursors to select the threshold value relative to the spike value');
            %[x_val_fund, y_val_fund] = ginput(1);
            % thresh_roi_fund = y_val_fund;
    
            %Use setdiff to find the difference. Isolates the spikes
            spikes_2 = setdiff(spikes, final, 'rows');
    
   
            %spikes_roi_2 = intersect(spikes_roi,final,'rows');
    
            %check to see if computation is correct
            %check = intersect(spikes_2, spikes_roi_2, 'rows'); %should return 0
    
            %count = 0;
            rf_filt_new2 = rf_filt_new(unwanted_row+1:end,:);
   
    
            %4) Remove spikes
            %replace all spikes occuring
            for i = 1:size(spikes_2,1)
                if(spikes_2(i,2)-1 == 0)
                    column_before = Col;
                else
                    column_before = spikes_2(i,2)-1;
                end
                if(spikes_2(i,2)+1 > Col)
                    column_after = 1;
                else
                    column_after = spikes_2(i,2) + 1;
                end
       
                vect_sum = rf_filt_new2(spikes_2(i,1), column_before) + rf_filt_new2(spikes_2(i,1),column_after);
       
                %vect_sum_fund = rf_filt_new2_fund(spikes_2_fund(i,1), column_before) + rf_filt_new2_fund(spikes_2_fund(i,1), column_after);
       
                vect_avg = vect_sum/2;
                %vect_avg_fund = vect_sum_fund/2;
       
                %replace with average
                rf_filt_new2(spikes_2(i,1), spikes_2(i,2)) = vect_avg;
                %rf_filt_new2_fund(spikes_2_fund(i,1), spikes_2_fund(i,2)) = vect_avg_fund;
                %count = count + 1;
            end
    
   
            %%first choose the background region without the spikes and the plaque
            %%without the spikes
             rf_filtered_new = abs(hilbert(rf_filt_new2));
            rf_filtered_new = rf_filtered_new/max(rf_filtered_new(:));
     
    
			
                %pad unwanted region with zeros
                %for i = 1:size(interest4,1)
                    %for j = 1:size(interest4,2)
                        %rf_filt_new(i,j) = 0;
                    %end
                %end
                %end of padding
            
             %figure;
             close(handle5);
             close(handle6);
            handle7 = figure('units','normalized','outerposition',[0 0 1 1])
            interest_new = [J K];
        
            T_new = 20*log10(rf_filtered_new);
      
            I_new = imagesc(20*log10(rf_filtered_new)); caxis([-60 -10]);
            colormap(gray);
            title('Select background for CTR computation. Be sure to select a region larger than or equal to the size of the ROI');
            pause
            BW_new = roipoly;
            [J_new, K_new] = find(BW_new > 0);
            interest_new_2 = [J_new K_new];
        
            for i = 1:length(J)
                interest_new(i,3) = T_new(interest_new(i,1), interest_new(i,2));
                interest_new(i,4) = T_new(interest_new_2(i,1), interest_new_2(i,2));
            end
	
            %get values of background region
            %for i = 1:length(J3)
            %interest_new(i,4) = T(interest3(i,1),interest3(i,2));
            %end
            %end of new code (2/10/14)
    
            %compute CTR
            plaque = mean(interest_new(:,3));
            background = mean(interest_new(:,4));
        
       
            CTR = plaque - background;
        
        
        	a = [plaque background CTR];
       
            %data(rows+counter,:) = a;%update data matrix
            %close all;
                %dlmwrite('data.txt',data)
                %clear all;
                %CTR
        end

        %recreate rf_filt_new
        rf_filt_new3 = zeros(unwanted_row,Col);
        rf_filt_new = [rf_filt_new3; rf_filt_new2];

        rf_filt_new = rf_filt_new(1:end,:);




        %generate final image
        rf_pad = zeros(start_pos_phase,Col);
        rf_filt_new_pad = [rf_pad' rf_filt_new']; 
        rf_filt_new_pad = rf_filt_new_pad./max(rf_filt_new_pad(:));

        figure; Polar2cart1(rf_filt_new_pad',Fs,1);
        caxis([-40 0]); title('Ultraharmonic image');
        disp('***************************************************');
        disp('FUNDAMENTAL IMAGE DATA');
        disp('***************************************************');
        disp(strcat('Fundamental CTR:',num2str(CTR_fund)));
        disp(strcat('Fundamental Background:',num2str(background_fund)));
        disp(strcat('Fundaemntal ROI:',num2str(plaque_fund)));
        disp('****************************************************');
        disp('');
        disp('****************************************************');
        disp('ULTRAHARMONIC IMAGE DATA');
        disp('*****************************************************');
        disp(strcat('Ultraharmonic CTR:',num2str(CTR)));
        disp(strcat('Ultraharmonic Background:',num2str(background)));
        disp(strcat('Ultraharmonic ROI:',num2str(plaque)));
        disp('*****************************************************');
        
    otherwise
        errordlg('Invalid selection. Exiting program','ERROR');
        error('Invalid selection');
       
end