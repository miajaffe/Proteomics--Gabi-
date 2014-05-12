raw_data = xlsread('Longitudinal_RawCounts_ForClass.xlsx', 2);

number_of_proteins = size(raw_data,1);
number_of_samples = size(raw_data,2);

raw_data_avg = mean(raw_data);

raw_data_avg_matrix = repmat(raw_data_avg, number_of_proteins, 1);

data_normalized = raw_data./raw_data_avg_matrix;

%% PCA 
% step 1: Calculate the average over all vectors
avg = mean(data_normalized);

%%
% step 2: Subtract the average from all vectors
diff_avg = data_normalized-repmat(avg, number_of_proteins, 1);

%%
% step 3: Calculate the covariance_matrixariance matrix

covariance_matrix = zeros(number_of_samples, number_of_samples);

for i = 1:number_of_samples;
    for j = 1:number_of_samples;
        covariance_matrix(i,j) = 1/(number_of_proteins - 1) * ...
            sum(diff_avg(:,i).*diff_avg(:,j));
    end
end
imagesc(covariance_matrix);
colorbar;
% disp(covariance_matrix);

%%
% step 4: calculate eigs of covariance_matrixariance matrix
[V,D] = eigs(covariance_matrix);
% disp(V);
% disp(diag(D));

%%
% plot data with the eigenvectors 
figure;
scatter3(data_normalized(:,1),data_normalized(:,2),...
    data_normalized(:,3),100,[1 0 0],'filled');
ev1 = V(:,1)*D(1,1);
ev2 = V(:,2)*D(2,2);
ev3 = V(:,3)*D(3,3);
hold on;
plot3([avg(1)-ev1(1) avg(1)+ev1(1)],[avg(2)-ev1(2) avg(2)+ev1(2)], ...
     [avg(3)-ev1(3) avg(3)+ev1(3)],'Color',[0 0 1],'LineWidth',3);
plot3([avg(1)-ev2(1) avg(1)+ev2(1)],[avg(2)-ev2(2) avg(2)+ev2(2)], ...
     [avg(3)-ev2(3) avg(3)+ev2(3)],'Color',[0 0 0],'LineWidth',3);
plot3([avg(1)-ev3(1) avg(1)+ev3(1)],[avg(2)-ev3(2) avg(2)+ev3(2)], ...
     [avg(3)-ev3(3) avg(3)+ev3(3)],'Color',[0 1 1],'LineWidth',3);
hold off;
axis equal

%%
 mapcaplot(data_normalized)






