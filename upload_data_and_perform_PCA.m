% This matlab code upload the vibrational data for six different colonies
% can calculates the PCA scores for the low brood and high brood normalised data.
% Written by Martin Bencsik, Sep. 2015.

clear

for averaging_factor = 30:30
    
    % Load the data for the channel 05:
    load(['C:\Users\SCI3BENCSM\Documents\Notes book\Swarmonitor\INRA_apiary_visits\matlab_code\new_apiary_2014_03_13\Publication_for_PS_alarm\brood_cycle_sensitivity\ch_05\without_normalisation\DFA_data'])
    % magenta:
    spectra_1 = data_brood_max;
    signal = mean(spectra_1);
    spectra_1 = spectra_1(:,find(signal < 2e-6));
    spectra_2 = data_brood_min;
    signal = mean(spectra_2);
    spectra_2 = spectra_2(:,find(signal < 1e-5));
    
    % Load the data for the channel 02:
    load(['C:\Users\SCI3BENCSM\Documents\Notes book\Swarmonitor\INRA_apiary_visits\matlab_code\new_apiary_2014_03_13\Publication_for_PS_alarm\brood_cycle_sensitivity\ch_02\without_normalisation\DFA_data'])
    
    temp_1 = data_brood_max;
    signal = mean(temp_1);
    spectra_1 = [spectra_1 temp_1(:,find(signal < 1.1e-6))];
    spectra_3 = data_brood_min;
    signal = mean(spectra_3);
    spectra_3 = spectra_3(:,find(signal < 4e-6));
    
    % Load the data for the channel 06:
    load(['C:\Users\SCI3BENCSM\Documents\Notes book\Swarmonitor\INRA_apiary_visits\matlab_code\new_apiary_2014_03_13\Publication_for_PS_alarm\brood_cycle_sensitivity\ch_06\without_normalisation\DFA_data'])
    
    temp_1 = data_brood_max;
    signal = mean(temp_1);
    spectra_1 = [spectra_1 temp_1(:,find(signal < 6e-7))];
    spectra_4 = data_brood_min;
    signal = mean(spectra_4);
    spectra_4 = spectra_4(:,find(signal < 2e-6));
    
    % Load the data for the channel 08:
    load(['C:\Users\SCI3BENCSM\Documents\Notes book\Swarmonitor\INRA_apiary_visits\matlab_code\new_apiary_2014_03_13\Publication_for_PS_alarm\brood_cycle_sensitivity\ch_08\without_normalisation\DFA_data'])
    
    temp_1 = data_brood_max;
    signal = mean(temp_1);
    spectra_1 = [spectra_1 temp_1(:,find(signal < 2e-6))];
    spectra_5 = data_brood_min;
    signal = mean(spectra_5);
    spectra_5 = spectra_5(:,find(signal < 1e-5));
    
    % Load the data for the channel 18:
    load(['C:\Users\SCI3BENCSM\Documents\Notes book\Swarmonitor\INRA_apiary_visits\matlab_code\new_apiary_2014_03_13\Publication_for_PS_alarm\brood_cycle_sensitivity\ch_18\without_normalisation\DFA_data'])
    
    temp_1 = data_brood_max;
    signal = mean(temp_1);
    spectra_1 = [spectra_1 temp_1(:,find(signal < 2e-6))];
    spectra_6 = data_brood_min;
    signal = mean(spectra_6);
    spectra_6 = spectra_6(:,find(signal < 1e-5));
    
    
    spectra_1 = spectra_1';
    spectra_2 = spectra_2';
    spectra_3 = spectra_3';
    spectra_4 = spectra_4';
    spectra_5 = spectra_5';
    spectra_6 = spectra_6';
    
    
    if averaging_factor > 1
        for sp_No = 1:6
            eval(['spectra_',num2str(sp_No),' = spectra_',num2str(sp_No),'(1:(averaging_factor*(floor(size(spectra_',num2str(sp_No),',1)/averaging_factor))),:);']);
            eval(['spectra_',num2str(sp_No),' = reshape(spectra_',num2str(sp_No),',averaging_factor,floor(size(spectra_',num2str(sp_No),',1)/averaging_factor),size(spectra_',num2str(sp_No),',2));']);
            eval(['spectra_',num2str(sp_No),' = squeeze(mean(spectra_',num2str(sp_No),'));']);
        end
    end
    
    
    data_set = [spectra_1; spectra_2; spectra_3; spectra_4; spectra_5; spectra_6];
    for uu = 1:size(data_set,1)
        data_set(uu,:) = (1/max(data_set(uu,:)))*data_set(uu,:);
    end
    for uu = 1:size(spectra_1,1)
        spectra_1(uu,:) = (1/max(spectra_1(uu,:)))*spectra_1(uu,:);
    end
    for uu = 1:size(spectra_2,1)
        spectra_2(uu,:) = (1/max(spectra_2(uu,:)))*spectra_2(uu,:);
    end
    for uu = 1:size(spectra_3,1)
        spectra_3(uu,:) = (1/max(spectra_3(uu,:)))*spectra_3(uu,:);
    end
    for uu = 1:size(spectra_4,1)
        spectra_4(uu,:) = (1/max(spectra_4(uu,:)))*spectra_4(uu,:);
    end
    for uu = 1:size(spectra_5,1)
        spectra_5(uu,:) = (1/max(spectra_5(uu,:)))*spectra_5(uu,:);
    end
    for uu = 1:size(spectra_6,1)
        spectra_6(uu,:) = (1/max(spectra_6(uu,:)))*spectra_6(uu,:);
    end
    
    
    % Centre the data set:
    temp_01 = mean(data_set);
    centred_data_set = (data_set - ones(size(data_set,1),1)*temp_01)';
    
    % Calculate the PCA scores and eigenvectors:
    L = centred_data_set*centred_data_set'; % L is the covariance matrix C=A*A'.
    [V D] = eig(L); % Diagonal elements of D are the eigenvalues for both L=A'*A and C=A*A'.
    
    centred_data_set = (spectra_1 - ones(size(spectra_1,1),1)*temp_01)';
    PCA_scores_01 = (V'*centred_data_set);
    PCA_scores_01 = flipud(PCA_scores_01);
    
    centred_data_set = (spectra_2 - ones(size(spectra_2,1),1)*temp_01)';
    PCA_scores_02 = (V'*centred_data_set);
    PCA_scores_02 = flipud(PCA_scores_02);
    
    centred_data_set = (spectra_3 - ones(size(spectra_3,1),1)*temp_01)';
    PCA_scores_03 = (V'*centred_data_set);
    PCA_scores_03 = flipud(PCA_scores_03);
    
    centred_data_set = (spectra_4 - ones(size(spectra_4,1),1)*temp_01)';
    PCA_scores_04 = (V'*centred_data_set);
    PCA_scores_04 = flipud(PCA_scores_04);
    
    centred_data_set = (spectra_5 - ones(size(spectra_5,1),1)*temp_01)';
    PCA_scores_05 = (V'*centred_data_set);
    PCA_scores_05 = flipud(PCA_scores_05);
    
    centred_data_set = (spectra_6 - ones(size(spectra_6,1),1)*temp_01)';
    PCA_scores_06 = (V'*centred_data_set);
    PCA_scores_06 = flipud(PCA_scores_06);
    
    Eigenspectra = fliplr(V);
    DFA_6D
      
end
%    print -dtiff -r300 Figure_07
