% clust.m
% Corban Swain, 2018

% clear the workspace
clear;

% load in data
data = load('cluster_data');

% annotate duplicate sample entries (with 'xx_n1', 'xx_n2', etc.)
nSampleNames = length(data.sample_names);
for iSamp = 1:nSampleNames
   name = data.sample_names{iSamp};
   dup_number = 0;
   for iSampTest = 1:nSampleNames
      if string(name) == string(data.sample_names{iSampTest})
         if dup_number == 0
            dup_number = dup_number + 1;
         elseif dup_number == 1
            data.sample_names{iSamp} = sprintf([name, '_n%d'], dup_number);
            dup_number = dup_number + 1;
            fprintf('Annotating a duplicate: %s \n', name)
         end
         
         if dup_number > 1
            data.sample_names{iSampTest} = sprintf([name, '_n%d'], ...
               dup_number);
            dup_number = dup_number + 1;
         end
      end
   end
end


%% Hierarchical Clustering
CNSUtils.FigureBuilder.setDefaults();
set(groot, 'defaultLineLineWidth', 0.5);

% use clustergram fuction to perform hierarchical clustering and 
% generate a dendrogram and heat map
cg = clustergram(data.data, 'RowLabels', data.sample_names, ...
   'ColumnLabels', data.gene_names, 'Colormap', redbluecmap);


%% k-Means Clustering
% setup constants
gene_k = 2;
sample_k = 6;
nReplicates = 200;

% perform k-means clustering along the sample axis
sample_idx = kmeans(data.data, sample_k, 'Replicates', nReplicates);
[ssi_values, sample_sort_idx] = sort(sample_idx);

% do the same along the gene axis
gene_idx = kmeans(data.data', gene_k, 'Replicates', nReplicates);
[gsi_values, gene_sort_idx] = sort(gene_idx);

% use the sorted indices to reorder the data matrices and label vectors
ordered_data = data.data(sample_sort_idx, :);
ordered_data = ordered_data(:, gene_sort_idx);
ord_gene_names = data.gene_names(gene_sort_idx);
ord_sample_names = data.sample_names(sample_sort_idx);

% plot the results
figure(1); clf;
h = heatmap(ord_gene_names, ord_sample_names, ordered_data);
h.ColorLimits = [-3, 3];
h.Colormap = redbluecmap;
h.GridVisible = 'off';
h.FontSize = 7;
fprintf('Done!\n');

% determine the endpoints of each cluster along the sample axis
n_sample = size(data.data, 1);
n_gene = size(data.data, 2);
sample_diffs = diff(ssi_values);
sample_branch_divisioins = ones(sample_k, 2);
sep_points = find(sample_diffs);
sample_branch_divisioins(1:(sample_k - 1), 2) = sep_points;
sample_branch_divisioins(2:sample_k, 1) = sep_points + 1;
sample_branch_divisioins(end) = n_sample

% and the gene axis
gene_diffs = diff(gsi_values);
gene_branch_divisioins = ones(gene_k, 2);
sep_points = find(gene_diffs);
gene_branch_divisioins(1:(gene_k - 1), 2) = sep_points;
gene_branch_divisioins(2:gene_k, 1) = sep_points + 1;
gene_branch_divisioins(end) = n_gene
