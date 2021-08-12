%%%%%%%%%%%%%% Graph Theory Measures %%%%%%%%%%%%
%%%%%%% Author: Rajanikant Panda
%%%%%%% Date of Development: 1st May 2017
%%%%%%% Date of Modification: 10th August 2021
%%%%%%% Supervised: Steven Laureys and Jitka Annen
%%%%%%% Reference papers:
%%%%%%% 1. Holla and Panda et al. (2017). Disrupted resting brain graph measures in individuals at high risk for alcoholism. Psychiatry Research: Neuroimaging, 265, 54-64.
%%%%%%% 2. Thibaut and Panda et al. (2021). Preservation of brain activity in unresponsive patients identifies MCS star. Annals of Neurology, 90(1), 89-100.
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% This Program compute the fallowing graph theory measures (absolute/real FC matrix) %%%%%%%%%%%%
%%%% Input = NxM matrix (N=No of ROIs, M=No of time point), this is the BOLD time series
%%%% Output = ROI to ROI connectivity matrix, Degree, clustering coeff (Network Segrigation)
%%%%            Participation coefficient (Network Integration), local and global efficiency
%%%%            path length, small-worldness
%%%% Reference Paper: 
%% %%%%%
clc
clear
%%% Provid the data folder path 
path='F:\CSG\fMRI\Meditation\GraphTheory\ToSent\GT_Post_Processing\Demo_Data\Control\';
cd(path)
SUBJlist=dir('CNT_*'); % write prefix of the group name
%%
for i=1:length(SUBJlist)
    SUBJname=SUBJlist(i).name;
    path1=([path SUBJname])
    cd(path1);
    %% 
    filelist= ([ SUBJname ]);   
   
        filename=filelist; 
        %prefix_name=filename(1:end);
                
        final_data=load([filename]);  % Loading the fMRI time series data
        final_data=final_data.y_roi_regressed; % Check that correct data is loadig %final_data.y_roi_regressed_filtered;

        GT_corr_data=corr(final_data);          % ROI to ROI connectivity matrix using pearsion correlations
        GT_corr_data_abs = abs(GT_corr_data);   % Making absolute to the connectivity matrix
        chanlocs = size(final_data,2);          % No  of ROIs      
        %% Sparcity based Thresholding
         sparsity_val=0.1:0.025:0.52; %sparsity_val=0.01:0.025:10.1:0.025:0.52; %sparsity_val=0.01:0.05:0.52; %sparsity_val=0.01:0.025:1;
        for i2=1:length (sparsity_val)
            %% %%%Network properties/ network measurement  for every sparsity threshold
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % calcualte binary matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf            
            %%
            GT_corr_data_thr(i2,:,:)=corr_data_thr; % asign to different array           
            GT_degree(i2,:)=degrees_und(corr_data_thr); % calculate degree             
            GT_clust_coeff(i2,:)=clustering_coef_bu(corr_data_thr); % clustering coeff (Network Segrigation)
            GT_local_eff(i2,:)=efficiency_bin(corr_data_thr,1); % local efficiency
            GT_global_eff(i2,:)=efficiency_bin(corr_data_thr,0); % global efficiency
            GT_distance_matrix(i2,:,:)=distance_bin(corr_data_thr); % distance matrix
            GT_path_length(i2)=charpath(squeeze(GT_distance_matrix(i2,:,:)),1,0); % path length
            %%%%% Participation coefficient (Network Integration)
            param.heuristic=50;            
            for i = 1:param.heuristic
                [Ci, allQ(i2,i)] = community_louvain(corr_data_thr);            
                allCi(i2,i,:) = Ci;  
                allpc(i2,i,:) = participation_coef(corr_data_thr,Ci); 
            end
            GT_modularity(i2)= mean(allQ(i2,:));  % modularity (Network Segrigation)
            GT_community_structure(i2,1:chanlocs) = squeeze(allCi(i2,1,:)); % community structure
            GT_participation_coeff(i2,1:chanlocs) = mean(squeeze(allpc(i2,:,:)));  %participation coefficient (Network Integration)
            %%%%%            
        end       
        varname=([SUBJname '_ABS'])
        save(varname);
    cd ..
end