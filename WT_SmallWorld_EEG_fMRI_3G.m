%% clear everything
clc; clear; close all

%% Group/Condition one

path = 'C:\Research\fMRI\GT\Group_1'
cd(path);
SUBJlist_Group1 = dir('Sub*');


%% Absolute and Random and Normalised CC, PL, SW for Pre_ET data extraction
for i = 1:length(SUBJlist_Group1)
    %%
    SUBJname = SUBJlist_Group1(i).name;
    path1=([path SUBJname]);
    cd(path1);
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group1_CC(i,:,:)= data.GT_clust_coeff;         
    Group1_PL(i,:,:)= data.GT_path_length;
    Group1_LE(i,:,:)= data.GT_local_eff;         
    Group1_GE(i,:,:)= data.GT_global_eff;
    Group1_Degree(i,:,:)= data.GT_degree;
    Group1_PC(i,:,:)= data.GT_participation_coeff;  
    
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group1_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group1_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group1_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group1_GE_rand(i,:,:,:)= data_rand.GT_global_eff_rand ;  %GT_global_eff_rand;
    Group1_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group1_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
end
%%
SpaRng=size(Group1_CC,2);
Group1_CC_50=Group1_CC(:,1:SpaRng,:);                          
Group1_CC_rand_squ = squeeze(mean(Group1_CC_rand,3));        
AvgGroup1_CC=mean(mean(Group1_CC_50,3));

sparsity_CC_Group1_50 = (mean(Group1_CC_50,3));
sparsity_CC_rand_Group1 = mean(Group1_CC_rand_squ,3);
sparsity_CC_rand_Group1_50 = sparsity_CC_rand_Group1 (:,1:SpaRng);

Group1_PL_50=Group1_PL(:,:,1:SpaRng);
AvgGroup1_PL=squeeze(mean(Group1_PL_50));

Sparsity_PL_Group1_50=squeeze(Group1_PL_50);
Sparsity_PL_Group1_rand_squ = squeeze(mean(Group1_PL_rand,3)); 
Sparsity_PL_Group1_rand_50 = Sparsity_PL_Group1_rand_squ (:,1:SpaRng);
Group1_Degree_50=Group1_Degree(:,1:SpaRng,:);                             
AvgGroup1_Degree=mean(mean(Group1_Degree_50,3));
sparsity_Degree_Group1_50 = (mean(Group1_Degree_50,3));

Group1_PC_50=Group1_PC(:,1:SpaRng,:);                          
Group1_PC_rand_squ = squeeze(mean(Group1_PC_rand,3));        
AvgGroup1_PC=mean(mean(Group1_PC_50,3));

sparsity_PC_Group1_50 = (mean(Group1_PC_50,3));
sparsity_PC_rand_Group1 = mean(Group1_PC_rand_squ,3);
sparsity_PC_rand_Group1_50 = sparsity_PC_rand_Group1 (:,1:SpaRng);

%%
for i = 1:length(SUBJlist_Group1); 
    for j = 1:SpaRng; 
            sparsity_CC_normalised_Group1(i,j) = sparsity_CC_Group1_50(i,j)/sparsity_CC_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = 1:SpaRng; 
            sparsity_PL_normalised_Group1(i,j) = Sparsity_PL_Group1_50(i,j)/Sparsity_PL_Group1_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:SpaRng; 
            SmallWorldNess_Group1(i,j) = sparsity_CC_normalised_Group1(i,j)/sparsity_PL_normalised_Group1(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:SpaRng; 
            sparsity_PC_normalised_Group1(i,j) = sparsity_PC_Group1_50(i,j)/sparsity_PC_rand_Group1_50(i,j); 
    end
end

%%
Group1_LE_50=Group1_LE(:,1:SpaRng,:);                          
Group1_LE_rand_squ = squeeze(mean(Group1_LE_rand,3));        
AvgGroup1_LE=mean(mean(Group1_LE_50,3));

sparsity_LE_Group1_50 = (mean(Group1_LE_50,3));
sparsity_LE_rand_Group1 = mean(Group1_LE_rand_squ,3);
sparsity_LE_rand_Group1_50 = sparsity_LE_rand_Group1 (:,1:SpaRng);

Group1_GE_50=Group1_GE(:,1:SpaRng);
AvgGroup1_GE=squeeze(mean(Group1_GE_50));

Sparsity_GE_Group1_50=squeeze(Group1_GE_50);
Sparsity_GE_Group1_rand_squ = squeeze(mean(Group1_GE_rand,3)); 
Sparsity_GE_Group1_rand_50 = Sparsity_GE_Group1_rand_squ (:,1:SpaRng); 

%%
for i = 1:length(SUBJlist_Group1); 
    for j = 1:SpaRng; 
            sparsity_LE_normalised_Group1(i,j) = sparsity_LE_Group1_50(i,j)/sparsity_LE_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = 1:SpaRng; 
            sparsity_GE_normalised_Group1(i,j) = Sparsity_GE_Group1_50(i,j)/Sparsity_GE_Group1_rand_50(i,j); 
    end
end




%% Group/Condition two study
path = 'C:\Research\fMRI\GT\Group_2'
cd(path);
SUBJlist_Group2 = dir('Sub*');
%%
for i = 1:length(SUBJlist_Group2)
    SUBJname = SUBJlist_Group2(i).name;
    path1=([path SUBJname]);
    cd(path1);
    %%
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group2_CC(i,:,:)= data.GT_clust_coeff;         
    Group2_PL(i,:,:)= data.GT_path_length;
    Group2_LE(i,:,:)= data.GT_local_eff;         
    Group2_GE(i,:,:)= data.GT_global_eff;
    Group2_Degree(i,:,:)= data.GT_degree;
    Group2_PC(i,:,:)= data.GT_participation_coeff;
    
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group2_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group2_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group2_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group2_GE_rand(i,:,:,:)= data_rand.GT_global_eff_rand;
    Group2_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group2_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
end
%%
Group2_CC_50=Group2_CC(:,1:SpaRng,:);                          
Group2_CC_rand_squ = squeeze(mean(Group2_CC_rand,3));        
AvgGroup2_CC=mean(mean(Group2_CC_50,3));

sparsity_CC_Group2_50 = (mean(Group2_CC_50,3));
sparsity_CC_rand_Group2 = mean(Group2_CC_rand_squ,3);
sparsity_CC_rand_Group2_50 = sparsity_CC_rand_Group2 (:,1:SpaRng);

Group2_PL_50=Group2_PL(:,:,1:SpaRng);
AvgGroup2_PL=squeeze(mean(Group2_PL_50));

Sparsity_PL_Group2_50=squeeze(Group2_PL_50);
Sparsity_PL_Group2_rand_squ = squeeze(mean(Group2_PL_rand,3)); 
Sparsity_PL_Group2_rand_50 = Sparsity_PL_Group2_rand_squ (:,1:SpaRng);
Group2_Degree_50=Group2_Degree(:,1:SpaRng,:);                             
AvgGroup2_Degree=mean(mean(Group2_Degree_50,3));
sparsity_Degree_Group2_50 = (mean(Group2_Degree_50,3));

Group2_PC_50=Group2_PC(:,1:SpaRng,:);                          
Group2_PC_rand_squ = squeeze(mean(Group2_PC_rand,3));        
AvgGroup2_PC=mean(mean(Group2_PC_50,3));

sparsity_PC_Group2_50 = (mean(Group2_PC_50,3));
sparsity_PC_rand_Group2 = mean(Group2_PC_rand_squ,3);
sparsity_PC_rand_Group2_50 = sparsity_PC_rand_Group2 (:,1:SpaRng);
%%
for i = 1:length(SUBJlist_Group2); 
    for j = 1:SpaRng; 
            sparsity_CC_normalised_Group2(i,j) = sparsity_CC_Group2_50(i,j)/sparsity_CC_rand_Group2_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group2); 
    for j = 1:SpaRng; 
            sparsity_PL_normalised_Group2(i,j) = Sparsity_PL_Group2_50(i,j)/Sparsity_PL_Group2_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:SpaRng; 
            SmallWorldNess_Group2(i,j) = sparsity_CC_normalised_Group2(i,j)/sparsity_PL_normalised_Group2(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:SpaRng; 
            sparsity_PC_normalised_Group2(i,j) = sparsity_PC_Group2_50(i,j)/sparsity_PC_rand_Group2_50(i,j); 
    end
end

%%
Group2_LE_50=Group2_LE(:,1:SpaRng,:);                          
Group2_LE_rand_squ = squeeze(mean(Group2_LE_rand,3));        
AvgGroup2_LE=mean(mean(Group2_LE_50,3));

sparsity_LE_Group2_50 = (mean(Group2_LE_50,3));
sparsity_LE_rand_Group2 = mean(Group2_LE_rand_squ,3);
sparsity_LE_rand_Group2_50 = sparsity_LE_rand_Group2 (:,1:SpaRng);

Group2_GE_50=Group2_GE(:,1:SpaRng);
AvgGroup2_GE=squeeze(mean(Group2_GE_50));

Sparsity_GE_Group2_50=squeeze(Group2_GE_50);
Sparsity_GE_Group2_rand_squ = squeeze(mean(Group2_GE_rand,3)); 
Sparsity_GE_Group2_rand_50 = Sparsity_GE_Group2_rand_squ (:,1:SpaRng); 

%%
for i = 1:length(SUBJlist_Group2); 
    for j = 1:SpaRng; 
            sparsity_LE_normalised_Group2(i,j) = sparsity_LE_Group2_50(i,j)/sparsity_LE_rand_Group2_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group2); 
    for j = 1:SpaRng; 
            sparsity_GE_normalised_Group2(i,j) = Sparsity_GE_Group2_50(i,j)/Sparsity_GE_Group2_rand_50(i,j); 
    end
end
%%  %%%%group 3 study %%%%%%
path = 'C:\Research\fMRI\GT\Group_3'
cd(path);
SUBJlist_Group3 = dir('Sub*');
%%
for i = 1:length(SUBJlist_Group3)
    SUBJname = SUBJlist_Group3(i).name;
    path1=([path SUBJname]);
    cd(path1);
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group3_CC(i,:,:)= data.GT_clust_coeff;         
    Group3_PL(i,:,:)= data.GT_path_length;
    Group3_LE(i,:,:)= data.GT_local_eff;         
    Group3_GE(i,:,:)= data.GT_global_eff;
    Group3_Degree(i,:,:)= data.GT_degree;
    Group3_PC(i,:,:)= data.GT_participation_coeff;
    
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group3_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group3_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group3_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group3_GE_rand(i,:,:,:)= data_rand.GT_global_eff_rand;
    Group3_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group3_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
end
%%
Group3_CC_50=Group3_CC(:,1:SpaRng,:);                          
Group3_CC_rand_squ = squeeze(mean(Group3_CC_rand,3));        
AvgGroup3_CC=mean(mean(Group3_CC_50,3));

sparsity_CC_Group3_50 = (mean(Group3_CC_50,3));
sparsity_CC_rand_Group3 = mean(Group3_CC_rand_squ,3);
sparsity_CC_rand_Group3_50 = sparsity_CC_rand_Group3 (:,1:SpaRng);

Group3_PL_50=Group3_PL(:,:,1:SpaRng);
AvgGroup3_PL=squeeze(mean(Group3_PL_50));

Sparsity_PL_Group3_50=squeeze(Group3_PL_50);
Sparsity_PL_Group3_rand_squ = squeeze(mean(Group3_PL_rand,3)); 
Sparsity_PL_Group3_rand_50 = Sparsity_PL_Group3_rand_squ (:,1:SpaRng);
Group3_Degree_50=Group3_Degree(:,1:SpaRng,:);                             
AvgGroup3_Degree=mean(mean(Group3_Degree_50,3));
sparsity_Degree_Group3_50 = (mean(Group3_Degree_50,3));

Group3_PC_50=Group3_PC(:,1:SpaRng,:);                          
Group3_PC_rand_squ = squeeze(mean(Group3_PC_rand,3));        
AvgGroup3_PC=mean(mean(Group3_PC_50,3));

sparsity_PC_Group3_50 = (mean(Group3_PC_50,3));
sparsity_PC_rand_Group3 = mean(Group3_PC_rand_squ,3);
sparsity_PC_rand_Group3_50 = sparsity_PC_rand_Group3 (:,1:SpaRng);
%%
for i = 1:length(SUBJlist_Group3); 
    for j = 1:SpaRng; 
            sparsity_CC_normalised_Group3(i,j) = sparsity_CC_Group3_50(i,j)/sparsity_CC_rand_Group3_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group3); 
    for j = 1:SpaRng; 
            sparsity_PL_normalised_Group3(i,j) = Sparsity_PL_Group3_50(i,j)/Sparsity_PL_Group3_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:SpaRng; 
            SmallWorldNess_Group3(i,j) = sparsity_CC_normalised_Group3(i,j)/sparsity_PL_normalised_Group3(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:SpaRng; 
            sparsity_PC_normalised_Group3(i,j) = sparsity_PC_Group3_50(i,j)/sparsity_PC_rand_Group3_50(i,j); 
    end
end

%%
Group3_LE_50=Group3_LE(:,1:SpaRng,:);                          
Group3_LE_rand_squ = squeeze(mean(Group3_LE_rand,3));        
AvgGroup3_LE=mean(mean(Group3_LE_50,3));

sparsity_LE_Group3_50 = (mean(Group3_LE_50,3));
sparsity_LE_rand_Group3 = mean(Group3_LE_rand_squ,3);
sparsity_LE_rand_Group3_50 = sparsity_LE_rand_Group3 (:,1:SpaRng);

Group3_GE_50=Group3_GE(:,1:SpaRng);
AvgGroup3_GE=squeeze(mean(Group3_GE_50));

Sparsity_GE_Group3_50=squeeze(Group3_GE_50);
Sparsity_GE_Group3_rand_squ = squeeze(mean(Group3_GE_rand,3)); 
Sparsity_GE_Group3_rand_50 = Sparsity_GE_Group3_rand_squ (:,1:SpaRng); 

%%
for i = 1:length(SUBJlist_Group3); 
    for j = 1:SpaRng; 
            sparsity_LE_normalised_Group3(i,j) = sparsity_LE_Group3_50(i,j)/sparsity_LE_rand_Group3_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group3); 
    for j = 1:SpaRng; 
            sparsity_GE_normalised_Group3(i,j) = Sparsity_GE_Group3_50(i,j)/Sparsity_GE_Group3_rand_50(i,j); 
    end
end

%% Absolute CC and PC ploting
figure (1)
y1_CC = mean (sparsity_CC_Group1_50);
z1_CC = std (sparsity_CC_Group1_50)/sqrt (length (sparsity_CC_Group1_50)); 
errorbar (y1_CC,z1_CC, 'k'); grid on; 
hold on
y2_CC = mean (sparsity_CC_Group2_50);
z2_CC = std (sparsity_CC_Group2_50)/sqrt (length (sparsity_CC_Group2_50));
errorbar (y2_CC,z2_CC, 'b'); grid on; 
y3_CC = mean (sparsity_CC_Group3_50);
z3_CC = std (sparsity_CC_Group3_50)/sqrt (length (sparsity_CC_Group3_50));
errorbar (y3_CC,z3_CC, 'r'); grid on; 
title('Absolute Clustering Coeficent')

figure (2)
y1_PC = mean (sparsity_PC_Group1_50);
z1_PC = std (sparsity_PC_Group1_50)/sqrt (length (sparsity_PC_Group1_50)); 
errorbar (y1_PC,z1_PC, 'k'); grid on; 
hold on
y2_PC = mean (sparsity_PC_Group2_50);
z2_PC = std (sparsity_PC_Group2_50)/sqrt (length (sparsity_PC_Group2_50));
errorbar (y2_PC,z2_PC, 'b'); grid on; 
y3_PC = mean (sparsity_PC_Group3_50);
z3_PC = std (sparsity_PC_Group3_50)/sqrt (length (sparsity_PC_Group3_50));
errorbar (y3_PC,z3_PC, 'r'); grid on; 
title('Absolute Participation Coeficent')
%%
figure (3)
y1_GE = mean (Sparsity_GE_Group1_50);
z1_GE = std (Sparsity_GE_Group1_50)/sqrt (length (Sparsity_GE_Group1_50)); 
errorbar (y1_GE,z1_GE, 'k'); grid on; 
hold on
y2_GE = mean (Sparsity_GE_Group2_50);
z2_GE = std (Sparsity_GE_Group2_50)/sqrt (length (Sparsity_GE_Group2_50));
errorbar (y2_GE,z2_GE, 'b'); grid on; 
y3_GE = mean (Sparsity_GE_Group3_50);

z3_GE = std (Sparsity_GE_Group3_50)/sqrt (length (Sparsity_GE_Group3_50));
errorbar (y3_GE,z3_GE, 'r'); grid on; 
title('Absolute Global Efficency')
%% %Ploting Normalized CC, PC, PL, GE, LE Images
figure (6)
y1_CC = mean (sparsity_CC_normalised_Group1);
z1_CC = std (sparsity_CC_normalised_Group1)/sqrt (length (sparsity_CC_normalised_Group1)); 
errorbar (y1_CC,z1_CC, 'k'); grid on; 
hold on
y2_CC = mean (sparsity_CC_normalised_Group2);
z2_CC = std (sparsity_CC_normalised_Group2)/sqrt (length (sparsity_CC_normalised_Group2));
errorbar (y2_CC,z2_CC, 'b'); grid on; 
y3_CC = mean (sparsity_CC_normalised_Group3);
z3_CC = std (sparsity_CC_normalised_Group3)/sqrt (length (sparsity_CC_normalised_Group3));
errorbar (y3_CC,z3_CC, 'r'); grid on; 
title('Normalised Clustering Coeficent')

figure (7)
hold on
y1_PL = mean (sparsity_PL_normalised_Group1);
z1_PL = std (sparsity_PL_normalised_Group1)/sqrt (length (sparsity_PL_normalised_Group1)); 
errorbar (y1_PL,z1_PL, 'k'); grid on; 
hold on
y2_PL = mean (sparsity_PL_normalised_Group2);
z2_PL = std (sparsity_PL_normalised_Group2)/sqrt (length (sparsity_PL_normalised_Group2));
errorbar (y2_PL,z2_PL, 'b'); grid on;
y3_PL = mean (sparsity_PL_normalised_Group3);
z3_PL = std (sparsity_PL_normalised_Group3)/sqrt (length (sparsity_PL_normalised_Group3));
errorbar (y3_PL,z3_PL, 'r'); grid on;
title('Normalised Path Length')

figure (8)
y1_SW = mean (SmallWorldNess_Group1);
z1_SW = std (SmallWorldNess_Group1)/sqrt (length (SmallWorldNess_Group1));
errorbar (y1_SW,z1_SW, 'k'); grid on;
hold on
y2_SW = mean (SmallWorldNess_Group2);
z2_SW = std (SmallWorldNess_Group2)/sqrt (length (SmallWorldNess_Group2));
errorbar (y2_SW,z2_SW, 'b'); grid on;
y3_SW = mean (SmallWorldNess_Group3);
z3_SW = std (SmallWorldNess_Group3)/sqrt (length (SmallWorldNess_Group3));
errorbar (y3_SW,z3_SW, 'r'); grid on;
title('Small Worldness')

figure (9)
y1_CC = mean (sparsity_LE_normalised_Group1);
z1_CC = std (sparsity_LE_normalised_Group1)/sqrt (length (sparsity_LE_normalised_Group1)); 
errorbar (y1_CC,z1_CC, 'k'); grid on; 
hold on
y2_CC = mean (sparsity_LE_normalised_Group2);
z2_CC = std (sparsity_LE_normalised_Group2)/sqrt (length (sparsity_LE_normalised_Group2));
errorbar (y2_CC,z2_CC, 'b'); grid on; 
y3_CC = mean (sparsity_LE_normalised_Group3);
z3_CC = std (sparsity_LE_normalised_Group3)/sqrt (length (sparsity_LE_normalised_Group3));
errorbar (y3_CC,z3_CC, 'r'); grid on; 
title('Normalised Local Eficency')

figure (10)
hold on
y1_PL = mean (sparsity_GE_normalised_Group1);
z1_PL = std (sparsity_GE_normalised_Group1)/sqrt (length (sparsity_GE_normalised_Group1)); 
errorbar (y1_PL,z1_PL, 'k'); grid on; 
hold on
y2_PL = mean (sparsity_GE_normalised_Group2);
z2_PL = std (sparsity_GE_normalised_Group2)/sqrt (length (sparsity_GE_normalised_Group2));
errorbar (y2_PL,z2_PL, 'b'); grid on;
y3_PL = mean (sparsity_GE_normalised_Group3);
z3_PL = std (sparsity_GE_normalised_Group3)/sqrt (length (sparsity_GE_normalised_Group3));
errorbar (y3_PL,z3_PL, 'r'); grid on;
title('Normalised Global Eficency')
figure (11)
y1_PC = mean (sparsity_PC_normalised_Group1);
z1_PC = std (sparsity_PC_normalised_Group1)/sqrt (length (sparsity_PC_normalised_Group1)); 
errorbar (y1_PC,z1_PC, 'k'); grid on;  
hold on
y2_PC = mean (sparsity_PC_normalised_Group2);
z2_PC = std (sparsity_PC_normalised_Group2)/sqrt (length (sparsity_PC_normalised_Group2));
errorbar (y2_PC,z2_PC, 'b'); grid on; 
y3_PC = mean (sparsity_PC_normalised_Group3);
z3_PC = std (sparsity_PC_normalised_Group3)/sqrt (length (sparsity_PC_normalised_Group3));
errorbar (y3_PC,z3_PC, 'r'); grid on; 
title('Normalised Participation Coeficent')
%%
Grp_con = ones(length(SUBJlist_Group1),1); Grp_mcs = 2*(ones(length(SUBJlist_Group2),1)); Grp_uws = 3*(ones(length(SUBJlist_Group3),1));
Grp = [Grp_con; Grp_mcs; Grp_uws]; 
Grp_PC = [trapz(sparsity_PC_normalised_Group1(:,10:10),2); trapz(sparsity_PC_normalised_Group2(:,10:10),2); trapz(sparsity_PC_normalised_Group3(:,10:10),2)];
figure(111); notBoxPlot(Grp_PC,Grp,0.5,'patch',ones(length(Grp_PC),1));
Grp_CC = [trapz(sparsity_CC_normalised_Group1(:,10:10),2); trapz(sparsity_CC_normalised_Group2(:,10:10),2); trapz(sparsity_CC_normalised_Group3(:,10:10),2)];
figure(112); notBoxPlot(Grp_CC,Grp,0.5,'patch',ones(length(Grp_CC),1));
Grp_GE = [trapz(sparsity_GE_normalised_Group1(:,10:10),2); trapz(sparsity_GE_normalised_Group2(:,10:10),2); mean(sparsity_GE_normalised_Group3(:,10:10),2)];
figure(113); notBoxPlot(Grp_GE,Grp,0.5,'patch',ones(length(Grp_GE),1));
%%
Grp_con = ones(length(SUBJlist_Group1),1); Grp_mcs = 2*(ones(length(SUBJlist_Group2),1)); Grp_uws = 3*(ones(length(SUBJlist_Group3),1));
Grp = [Grp_con; Grp_mcs; Grp_uws]; 
Grp_PC = [(sparsity_PC_normalised_Group1(:,10:10));(sparsity_PC_normalised_Group2(:,10:10)); (sparsity_PC_normalised_Group3(:,10:10))];
figure(111); notBoxPlot(Grp_PC,Grp,0.5,'patch',ones(length(Grp_PC),1));
Grp_CC = [(sparsity_CC_normalised_Group1(:,10:10)); (sparsity_CC_normalised_Group2(:,10:10)); (sparsity_CC_normalised_Group3(:,10:10))];
figure(112); notBoxPlot(Grp_CC,Grp,0.5,'patch',ones(length(Grp_CC),1));
Grp_GE = [(sparsity_GE_normalised_Group1(:,10:10)); (sparsity_GE_normalised_Group2(:,10:10)); mean(sparsity_GE_normalised_Group3(:,10:10))];
figure(113); notBoxPlot(Grp_GE,Grp,0.5,'patch',ones(length(Grp_GE),1));
%% --------------t-stats between Group2 and Group1 sparsity level------------%

for i = 1:SpaRng
    [h_CC1(i),p_CC1(i)] = ttest(sparsity_CC_normalised_Group2(:,i),sparsity_CC_normalised_Group1(:,i),0.05,'right');
end

h_CC1
p_CC1
%%
for i = 1:SpaRng
    [h_GE1(i),p_GE1(i)] = ttest2(sparsity_GE_normalised_Group3(:,i),sparsity_GE_normalised_Group2(:,i),0.05,'left');
end
h_GE1
p_GE1
%%
for i = 1:SpaRng
    [h_SW1(i),p_SW1(i)] = ttest2(SmallWorldNess_Group2(:,i),SmallWorldNess_Group1(:,i),0.05,'right');
end

h_SW1
p_SW1
%%
for i = 1:SpaRng
    [h_PC1(i),p_PC1(i)] = ttest(sparsity_PC_normalised_Group1(:,i),sparsity_PC_normalised_Group2(:,i),0.05,'left');
end
h_PC1
p_PC1
%%
%----------------------%% Brain resion significant Computations for CC %-------------%


sparsity_CC_Group1_ROI = squeeze(mean(Group1_CC_50,2)); 
sparsity_CC_Group2_ROI = squeeze(mean(Group2_CC_50,2));
sparsity_CC_Group3_ROI = squeeze(mean(Group3_CC_50,2));
sparsity_CC_rand_Group1_ROI = mean(Group1_CC_rand_squ,2);
sparsity_CC_rand_Group2_ROI = mean(Group2_CC_rand_squ,2);
sparsity_CC_rand_Group3_ROI = mean(Group3_CC_rand_squ,2);


for i = 1:length(SUBJlist_Group1); 
    for j = 1:268; 
            sparsity_CC_normalised_Group1_ROI(i,j) = sparsity_CC_Group1_ROI(i,j)/sparsity_CC_rand_Group1_ROI(i,j);
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:268; 
            sparsity_CC_normalised_Group2_ROI(i,j) = sparsity_CC_Group2_ROI(i,j)/sparsity_CC_rand_Group2_ROI(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:268; 
            sparsity_CC_normalised_Group3_ROI(i,j) = sparsity_CC_Group3_ROI(i,j)/sparsity_CC_rand_Group3_ROI(i,j); 
    end
end
%%
%----------------------%% Brain resion significant Computations for PC %-------------%

sparsity_PC_Group1_ROI = squeeze(mean(Group1_PC_50,2)); 
sparsity_PC_Group2_ROI = squeeze(mean(Group2_PC_50,2)); 
sparsity_PC_Group3_ROI = squeeze(mean(Group3_PC_50,2));
sparsity_PC_rand_Group1_ROI = mean(Group1_PC_rand_squ,2);
sparsity_PC_rand_Group2_ROI = mean(Group2_PC_rand_squ,2);
sparsity_PC_rand_Group3_ROI = mean(Group3_PC_rand_squ,2);



for i = 1:length(SUBJlist_Group1); 
    for j = 1:268; 
            sparsity_PC_normalised_Group1_ROI(i,j) = sparsity_PC_Group1_ROI(i,j)/sparsity_PC_rand_Group1_ROI(i,j);
    end
end



for i = 1:length(SUBJlist_Group2); 
    for j = 1:268; 
            sparsity_PC_normalised_Group2_ROI(i,j) = sparsity_PC_Group2_ROI(i,j)/sparsity_PC_rand_Group2_ROI(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:268; 
            sparsity_PC_normalised_Group3_ROI(i,j) = sparsity_PC_Group3_ROI(i,j)/sparsity_PC_rand_Group3_ROI(i,j); 
    end
end
%% Brain resion significant Computations for PC

for i = 1:268
    [h_CC_ROI_Normalised(i),p_CC_ROI_Normalised(i)] = ttest(sparsity_CC_normalised_Group1_ROI(:,i),sparsity_CC_normalised_Group2_ROI(:,i),0.01,'right');
end

h_CC_ROI_Normalised
p_CC_ROI_Normalised
%%
clc
for i = 1:268
    [h_PC_ROI_Normalised(i),p_PC_ROI_Normalised(i)] = ttest(sparsity_PC_normalised_Group1_ROI(:,i),sparsity_PC_normalised_Group2_ROI(:,i),0.01,'left');
end

h_PC_ROI_Normalised   
p_PC_ROI_Normalised
[p,h]=fdr(p_CC_ROI_Normalised,0.05);
p

%%  Graphical Plot od GT measures %%%
% % % tpz_cc_control = trapz(sparsity_CC_normalised_Control(:,1:30),2)/30;
% % % tpz_cc_Patient = trapz(sparsity_CC_normalised_Patient(:,1:30),2)/30;
% % % tpz_pl_control = trapz(sparsity_PL_normalised_Control(:,1:30),2)/30;
% % % tpz_pl_Patient = trapz(sparsity_PL_normalised_Patient(:,1:30),2)/30;
% % % tpz_sw_control = trapz(SmallWorldNess_Control(:,1:30),2)/30;
% % % tpz_sw_Patient = trapz(SmallWorldNess_Patient(:,1:30),2)/30;
% % % tpz_le_control = trapz(sparsity_LE_normalised_Control(:,1:30),2)/30;
% % % tpz_le_Patient = trapz(sparsity_LE_normalised_Patient(:,1:30),2)/30;
% % % tpz_ge_control = trapz(sparsity_GE_normalised_Control(:,1:30),2)/30;
% % % tpz_ge_Patient = trapz(sparsity_GE_normalised_Patient(:,1:30),2)/30;
% % % 
% % % tpz_GT = [tpz_cc_control tpz_cc_Patient tpz_pl_control tpz_pl_Patient tpz_sw_control tpz_sw_Patient tpz_le_control tpz_le_Patient tpz_ge_control tpz_ge_Patient]
% % % 
% % % [p,h]=ttest2(tpz_GT(:,1),tpz_GT(:,2), 0.05,'left')
% % % 
% % % 
