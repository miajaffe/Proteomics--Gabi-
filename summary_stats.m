%%Summary Stats

clear all
close all hidden
clc
load('axes140523.mat');
load('LetterMap.mat');
load('normOverlordFinal_140523.mat');
load('ProteinMap.mat');
load('MusProt.mat')
location_labels = axes{4};
colonization_labels = axes{3};
normOverlord = normOverlordFinal * 100; 
proteins = axes{1}'
% for location = 1:5
%         replicate_sum = sum(normOverlord(:,:,:,location), 2);
%         total_sum = sum(replicate_sum, 3);
%         total_prot_abundance = sum(total_sum, 1);
%         n_proteins = length(total_sum);
%         total_abundance_vect = repmat(total_prot_abundance, n_proteins, 1);
%         prot_percentages = total_sum ./ total_abundance_vect;
%         [rows, cols] = find(prot_percentages < 0.01);
%         prot_small_abundance = total_sum(rows,1);
%         summed_small_abundance = sum(prot_small_abundance);
%         [rows2, cols2, vals2] = find(prot_percentages > 0.01);
%         prot_large_abundance = total_sum(rows2,1);
%         final_abundances = [prot_large_abundance;summed_small_abundance]; 
%         subplot(2,3,location)
%         pie(final_abundances)
%         title(location_labels(location));
% end

%same thing but separate by location 
 final_table = {};
 final_names = {};
 count = 1;      
 cutoff = 0.02; %cutoff protein abundance to group into other
for location = 1:5
    for colonization = 1:3
        title_label = strcat(location_labels(location), colonization_labels(colonization));
        % sum up protein abundance for all of the replicates
        replicate_sum = sum(normOverlord(:,:,colonization,location), 2); 
        % total protein abundance for all proteins 
        total_prot_abundance = sum(replicate_sum, 1);
        n_proteins = length(replicate_sum);
        total_abundance_vect = repmat(total_prot_abundance, n_proteins, 1);
        %calculate the percentage of each protein abundance compared to the
        %total
        prot_percentages = replicate_sum ./ total_abundance_vect;
        % Group all of the proteins with < 1% abundance into one other
        % group
        [rows, cols] = find(prot_percentages <= cutoff);
        prot_small_abundance = replicate_sum(rows,1);
        summed_small_abundance = sum(prot_small_abundance);
        [rows2, cols2] = find(prot_percentages > cutoff);
        fprintf('%d', count);
        large_prot = proteins(rows2,1)
        prot_large_abundance = replicate_sum(rows2,1)
        prot_name_counts = {large_prot prot_large_abundance}
        final_table = [final_table prot_name_counts];
        final_abundances = [prot_large_abundance; summed_small_abundance]; 
        figure 
        %%subplot(5,3,count)
        abundance = pie(final_abundances)
        title(title_label(1));
        count = count + 1;
        % saveas(abundance,strcat('figure', num2str(count), '.pdf'))
        
        
        %set colors
        hp = findobj(abundance, 'Type', 'patch');
        for i = 1: length(large_prot)
            switch(large_prot{i})
                case 'A8DUK4'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'B2RS76'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'P02815'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'P05208'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'P07146'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'P07743'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q3SYP2'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q61900'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q62395'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q64097'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q6P8U6'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q7TNM8'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q7TPZ8'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q8C5B4'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q8K0C5'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q8R0I0'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q91WB5'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q91WL7'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q91X79'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case'Q91XA9'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q921Y7'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9CPN9'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9CQ52'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9CQC2'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9CR35'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9CY06'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9ER05'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9JM84'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9R0T7'
                    set(hp(i), 'FaceColor', [0 1 1]);
                case 'Q9CQC2'  
                    set(hp(i), 'FaceColor', [0 1 1]);
                   
                    
            end
        end
    end
end


