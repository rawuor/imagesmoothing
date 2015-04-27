function [spikes, rf_filt_new3, rf_filt_new2] = no_ROI_ultraharmonic(rf_filt_new,start_pos_phase, Fs, Col)
    rf_filtered = abs(hilbert(rf_filt_new));
    rf_filtered = rf_filtered/max(rf_filtered(:));

    handle1 = figure('units','normalized','outerposition', [0 0 1 1]);
    imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
    colormap(gray);
    title('Select the unwanted region');
    pause
    BW4 = roipoly;
    [J4, K4] = find(BW4 > 0); %find the coordinates of the unwanted region
    interest4 = [J4 K4];
    
    %Remove unwanted region
    unwanted_row = max(interest4(:,1));
    rf_filtered = rf_filtered(unwanted_row+1:end, :);
    
    close(handle1); %not needed
    handle2 = figure('units','normalized','outerposition', [0 0 1 1]);
    imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
    colormap(gray);
    title('Use cursors to select any spike.');
    [spike_col, spike_row] = ginput(1); 
    
    title('Now select the background region. Be sure not to select a spiked region');
    pause
    BW3 = roipoly;
    [J3 K3] = find(BW3 > 0); %find the coordinates of the background region
    interest3 = [J3 K3]; %background region
    
    %get values fo the background region
    for i = 1:size(interest3,1)
        for j = 1:size(interest3, 2)
            interest3(i,3) = rf_filtered(interest3(i,1), interest3(i,2));
        end
    end
    thresh = max(interest3(:,3)); %threshold value for the background region
    
    %evaluate the spike relative to the background region
    handle3 = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    plot(rf_filtered(:,round(spike_col)));
    string = num2str(round(spike_row));
    string = strcat('Use the cursors to select the threshold value relative to the spike value. Your spike is on x-value =', string);
    title(string);
    [x_val, y_val] = ginput(1);
    thresh_roi = y_val;
    
    [row, col] = find(rf_filtered > thresh);
    
    %roi threshold selection 
    [row_new, col_new] = find(rf_filtered > thresh_roi);
    spikes_roi = [row_new, col_new];
    
    spikes = [row col];
    
    spikes_2 = spikes;
    
    rf_filt_new2 = rf_filt_new(unwanted_row+1:end,:);
    rf_filt_new3 = rf_filt_new(unwanted_row+1:end,:);
    
    
    %Remove spikes
   for i = 1:size(spikes_2,1)
       if (spikes_2(i,2) - 1 == 0)
           column_before = Col;
       else
           column_before = spikes_2(i,2) - 1;
       end
       if(spikes_2(i,2) + 1 > Col)
           column_after = 1;
       else
           column_after = spikes_2(i,2) + 1;
       end
       
       vect_sum = rf_filt_new2(spikes_2(i,1), column_before) + rf_filt_new2(spikes_2(i,1), column_after);
       
       vect_avg = vect_sum/2;
       
       rf_filt_new2(spikes_2(i,1), spikes_2(i,2)) = vect_avg;
   end
    rf_filtered_new = abs(hilbert(rf_filt_new2));
    rf_filtered_interim = abs(hilbert(rf_filt_new3));
    close(handle2);
    close(handle3);
    rf_filtered_new = rf_filtered_new./max(rf_filtered_interim(:));
    handle4 = figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(20*log10(rf_filtered_new)); caxis([-60 -10]);
    colormap(gray);
    
    %Recreate rf_filt_new
    rf_filt_new4 = zeros(unwanted_row, Col);
    rf_filt_new5 = [rf_filt_new4; rf_filt_new2];
    
    rf_filt_new5 = rf_filt_new5(1:end, :);
    
    %generate final image
    rf_pad = zeros(start_pos_phase, Col);
    rf_filt_new_pad = [rf_pad' rf_filt_new5'];
    rf_filt_new_pad = rf_filt_new_pad./max(rf_filtered_interim(:));
    close(handle4);
    figure;
    Polar2cart1(rf_filt_new_pad', Fs, 1);
    caxis([-40 0]); title('Ultraharmonic image');
    
    %CTR computation
    background = 20*log10(mean(interest3(:,3)));
    plaque = background;
    CTR = plaque - background;
    disp('*************************************************');
    disp(strcat('CTR:',num2str(CTR)));
    disp(strcat('Background:',num2str(background)));
    disp(strcat('ROI (Your ROI has the same value as your background):',num2str(plaque)));
    disp('***************************************************');
    
end