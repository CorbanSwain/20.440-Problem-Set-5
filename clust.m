% clust.m
% Corban Swain, 2018
clear;
data = load('cluster_data');

% remove duplicates
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
cg = clustergram(data.data, 'RowLabels', data.sample_names, ...
   'ColumnLabels', data.gene_names, 'Colormap', redbluecmap);

%% k-Means Clustering
gene_k = 3;
sample_k = 5;
nReplicates = 200;
sample_idx = kmeans(data.data, 5, 'Replicates', nReplicates);
[~, sample_sort_idx] = sort(sample_idx);
gene_idx = kmeans(data.data', sample_k, 'Replicates', nReplicates);
[~, gene_sort_idx] = sort(gene_idx);

ordered_data = data.data(sample_sort_idx, :);
ordered_data = ordered_data(:, gene_sort_idx);
ord_gene_names = data.gene_names(gene_sort_idx);
ord_sample_names = data.sample_names(sample_sort_idx);

figure(1); clf;
h = heatmap(ord_gene_names, ord_sample_names, ordered_data);
h.ColorLimits = [-3, 3];
h.Colormap = redbluecmap;
h.GridVisible = 'off';
h.FontSize = 5;
fprintf('Done!\n');

