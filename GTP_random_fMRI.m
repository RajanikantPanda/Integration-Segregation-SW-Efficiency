%%%%%%%%%%%%%% Graph Theory Measures %%%%%%%%%%%%
%%%%%%% Author: Rajanikant Panda
%%%%%%% Date of Development: 1st May 2017
%%%%%%% Date of Modification: 10th August 2021
%%%%%%% Supervised: Steven Laureys and Jitka Annen
%%%%%%% Reference papers:
%%%%%%% 1. Holla and Panda et al. (2017). Disrupted resting brain graph measures in individuals at high risk for alcoholism. Psychiatry Research: Neuroimaging, 265, 54-64.
%%%%%%% 2. Thibautand Panda et al. (2021). Preservation of brain activity in unresponsive patients identifies MCS star. Annals of Neurology, 90(1), 89-100.
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% This Program compute the fallowing graph theory measures (for random FC matrix by conening same no of node and links as absolute/original FC) %%%%%%%%%%%%
%%%% Input = NxM matrix (N=No of ROIs, M=No of time point), this is the BOLD time series
%%%% Output =   Random Degree, Random clustering coeff (Network Segrigation)
%%%%            Random Participation coefficient (Network Integration), Random local and global efficiency
%%%%            Random path length, Random small-worldness

%% %%%%%
clc
clear
%%% Provid the data folder path 
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
%        prefix_name=filename(1:end);
                
        final_data=load([filename]);   % Loading the fMRI time series data
        final_data=final_data.y_roi_regressed; % Check that correct data is loadig %final_data.y_roi_regressed_filtered;

        GT_corr_data=corr(final_data);          % ROI to ROI connectivity matrix using pearsion correlations
        GT_corr_data_abs = abs(GT_corr_data);   % Making absolute to the connectivity matrix
        chanlocs = size(final_data,2);          % No  of ROIs                       
        %% Sparcity based Thresholding
        sparsity_val=0.1:0.025:0.52; %sparsity_val=0.01:0.05:0.52; %sparsity_val=0.01:0.025:1;
        for i2=1:length (sparsity_val)
            %% %%%Network properties/ network measurement for every sparsity threshold
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % calcualte binary matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf
                       
            for random_number=1:40   % No of Randomization Matrix to be generate
                             
                random_network=randmio_und(corr_data_thr,5);
                %%
                GT_corr_data_rand_thr(i2,random_number,:,:)=random_network;% asign to different array
                GT_degree_rand(i2,random_number,:)=degrees_und(random_network); % calculate degree
                GT_clust_coeff_rand(i2,random_number,:)=clustering_coef_bu(random_network);%% clustering coeff (Network Segrigation)            
                GT_local_eff_rand(i2,random_number,:)=efficiency_bin(random_network,1); % local efficiency 
                GT_global_eff_rand(i2,random_number,:)=efficiency_bin(random_network,0); % global efficiency
                GT_distance_matrix_rand(i2,random_number,:,:)=distance_bin(random_network); % distance matrix
                GT_path_length_rand(i2,random_number) =charpath(squeeze(GT_distance_matrix_rand(i2,random_number,:,:)),1,0); % path length
            
            %%%%% Participation coefficient (Network Integration)
                param.heuristic=50;
                for i = 1:param.heuristic
                    [Ci, allQ(i2,random_number,i)] = community_louvain(random_network);
                     allCi(i2,random_number,i,:) = Ci;
                     allpc(i2,random_number,i,:) = participation_coef(random_network,Ci); 
                end
                GT_modularity_rand(i2,random_number)= mean(allQ(i2,random_number,:));  % modularity
                GT_community_structure_rand(i2,random_number,1:chanlocs) = squeeze(allCi(i2,random_number,1,:)); % community structure
                GT_participation_coeff_rand(i2,random_number,1:chanlocs) = mean(squeeze(allpc(i2,random_number,:,:)));  %participation coefficient (Network Integration)
                %%%%%    
            end
        end
        varname=([SUBJname '_RAND'])
        save(varname);
    cd ..
end
