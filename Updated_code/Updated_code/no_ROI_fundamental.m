function no_ROI_fundamental(rf_filt_new, rf_neg_phase, start_pos_phase, Fs, Col)
  rf_filtered = abs(hilbert(rf_filt_new));
  rf_filtered = rf_filtered./max(rf_filtered(:));
  
  %Prepare fundamental matrix for spike removal
  rf_filt_new_fund = rf_neg_phase;
  
  %normalization
  rf_filtered_fund = abs(hilbert(rf_filt_new_fund));
  rf_filtered_fund = rf_filtered_fund/max(rf_filtered_fund(:));
  
  %display image to get rid of unwanted middle part
  handle1 = figure('units','normalized','outerposition',[0 0 1 1]);
  imagesc(20*log10(rf_filtered)); caxis([-60 -10]);
  colormap(gray);
  title('Select unwanted region');
  pause
  BW4 = roipoly;
  [J4, K4] = find(BW4 > 0); %Find the coordinates of the unwanted region
  interest4 = [J4, K4];
  
  %Remove the unwanted region
  unwanted_row = max(interest4(:,1));
  rf_filtered = rf_filtered(unwanted_row+1:end,:);
  rf_filtered_fund = rf_filtered_fund(unwanted_row+1:end,:);
  
  close(handle1);
  handle2 = figure('units','normalized','outerposition',[0 0 1 1]);
  imagesc(20*log10(rf_filtered_fund)); caxis([-60 4]);
  colormap(gray);
  title('Use cursors to select a spike.');
  [spike_col_fund, spike_row_fund] = ginput(1);
  
  title('Select the background region. Be sure not to select a spiked region');
  pause 
  BW3_fund = roipoly;
  [J3_fund, K3_fund] = find(BW3_fund > 0);
  interest3_fund = [J3_fund, K3_fund]; %region of interest
  
  %get values of background region
  for i = 1:size(interest3_fund,1)
    for j = 1:size(interest3_fund,2)
        interest3_fund(i,3) = rf_filtered_fund(interest3_fund(i,1), interest3_fund(i,2));
    end
  end
  thresh_fund = max(interest3_fund(:,3)); %threshold value of the background region
  
  %evaluate the spike relative to the background region
  handle3 = figure('units','normalized','outerposition',[0 0 1 1]);
  hold on
  plot(rf_filtered_fund(:,round(spike_col_fund)));
  string = strcat('Use the cursors to select the threshold value relative to the spike value. Your spike value is on x-value = ',num2str(round(spike_row_fund)));
  title(string);
  
  [x_val_fund, y_val_fund] = ginput(1);
  thresh_roi_fund = y_val_fund;
  
  [row_fund, col_fund] = find(rf_filtered_fund > thresh_fund);
  
  %roi threshold selection
  [row_new_fund, col_new_fund] = find(rf_filtered_fund > thresh_roi_fund);
  spikes_roi_fund = [row_new_fund, col_new_fund];
  
  spikes_fund = [row_fund, col_fund];
  
  rf_filt_new2_fund = rf_filt_new_fund(unwanted_row+1:end,:);
  
  %Remove fundamental image spikes
  for j = 1:size(spikes_fund,1)
    if(spikes_fund(j,2)- 1 == 0)
        column_before_fund = Col;
    else
        column_before_fund = spikes_fund(j,2) - 1;
    end
    
    if(spikes_fund(j,2) + 1 > Col)
        column_after_fund = 1;
    else
        column_after_fund = spikes_fund(j,2) + 1;
    end
    
    vect_sum_fund = rf_filt_new2_fund(spikes_fund(j,1), column_before_fund) + rf_filt_new2_fund(spikes_fund(j,1), column_after_fund);
    vect_avg_fund = vect_sum_fund/2;
    rf_filt_new2_fund(spikes_fund(j,1), spikes_fund(j,2)) = vect_avg_fund; %assign averaged value
  end
  
  rf_filtered_new_fund = abs(hilbert(rf_filt_new2_fund));
  rf_filtered_new_fund = rf_filtered_new_fund./max(rf_filtered_new_fund(:));
  %rf_filtered_new_fund = rf_filtered_new_fund/max(abs(hilbert(rf_filt_new_fund)));
  close(handle2);
  close(handle3);
  handle4 = figure('units','normalized','outerposition',[0 0 1 1]);
  imagesc(20*log10(rf_filtered_new_fund)); caxis([-60 4]);
  colormap(gray);
  
  background_fund = 20*log10(mean(interest3_fund(:,3)));
  plaque_fund = background_fund;
  CTR = plaque_fund - background_fund;
  
  disp('******************************************************************************************');
  disp(strcat('CTR:',num2str(CTR)));
  disp(strcat('Background:',num2str(background_fund)));
  disp(strcat('ROI (value is equal to background since there is no ROI):',num2str(plaque_fund)));
  disp('******************************************************************************************');
  %recreate rf_filt_new_fund
  rf_filt_new3_fund =zeros(unwanted_row,Col);
 
  rf_filt_new_fund = [rf_filt_new3_fund; rf_filt_new2_fund];
  
  rf_filt_new_fund = rf_filt_new_fund(1:end,:);
  
  %generate final fundamental iamge
  rf_pad_fund = zeros(start_pos_phase, Col);
  rf_filt_new_pad_fund = [rf_pad_fund; rf_filt_new_fund];
  rf_filt_new_pad_fund = rf_filt_new_pad_fund./max(rf_filt_new_pad_fund(:));
  
  size(rf_filt_new_pad_fund)
  close(handle4);
  figure; Polar2cart1(rf_filt_new_pad_fund, Fs, 1);
  caxis([-60 4])
  title('Fundamental image');
  
 
end