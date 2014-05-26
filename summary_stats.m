%%Summary Stats

clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('normOverlordF.mat');
load('ProteinMap.mat');
load('MusProt.mat')
location_labels = axes{4};
colonization_labels = axes{3};
normOverlord = normOverlordF * 100; 
for location = 1:5
        replicate_sum = sum(normOverlord(:,:,:,location), 2);
        total_sum = sum(replicate_sum, 3);
        total_prot_abundance = sum(total_sum, 1);
        n_proteins = length(total_sum);
        total_abundance_vect = repmat(total_prot_abundance, n_proteins, 1);
        prot_percentages = total_sum ./ total_abundance_vect;
        [rows, cols] = find(prot_percentages < 0.01);
        prot_small_abundance = total_sum(rows,1);
        summed_small_abundance = sum(prot_small_abundance);
        [rows2, cols2, vals2] = find(prot_percentages > 0.01);
        prot_large_abundance = total_sum(rows2,1);
        final_abundances = [prot_large_abundance;summed_small_abundance]; 
        subplot(2,3,location)
        pie(final_abundances)
        title(location_labels(location));
end

%same thing but separate by location 
 count = 1;       
for location = 1:5
    for colonization = 1:3
        replicate_sum = sum(normOverlord(:,:,colonization,location), 2);
        total_prot_abundance = sum(replicate_sum, 1);
        n_proteins = length(replicate_sum);
        total_abundance_vect = repmat(total_prot_abundance, n_proteins, 1);
        prot_percentages = replicate_sum ./ total_abundance_vect;
        [rows, cols] = find(prot_percentages < 0.01);
        prot_small_abundance = total_sum(rows,1);
        summed_small_abundance = sum(prot_small_abundance);
        [rows2, cols2, vals2] = find(prot_percentages > 0.01);
        prot_large_abundance = total_sum(rows2,1);
        final_abundances = [prot_large_abundance;summed_small_abundance]; 
        %subplot(5,3,count)
        figure
        pie(final_abundances)
        title_label = strcat(location_labels(location), colonization_labels(colonization));
       
        title(title_label(1));
        count = count + 1;
    end
end

