%% settings
c=flipud(gray);
set(0,'DefaultFigureColormap',c)

%% read data from csv, convert to struct
T = readtable("data_reduced.csv","ReadRowNames",true); % don't type T in the command window. The display is super slow!!
metadata = readtable("metadata_reduced.csv","ReadRowNames",true);
T_merged=join(T,metadata,'LeftKeys','Unnamed_0','RightKeys','x1_ID'); % don't type T_merged in the command window. The display is super slow!!
data=struct;
data.cellNames = T.Properties.RowNames;
data.geneNames = T.Properties.VariableNames(2:end);
data.clust = T_merged{:,'x7_clust_all'};
data.data = T{1:end,1:end};
clear T T_merged


%% 1. normalize data
% Now you have the raw data in data.data, which is a n_cell x n_gene matrix.
% First normalize the data such that the analysis won't be dominated by a
% few genes with very high expression counts. 
% There are many different ways to normalize, but we will normalize such
% that the standard deviation of each columns is 1.
% We will work with the normalized data for all the later sections.



%% 2. visualize data
% Visualize the normalized data. You might use imagesc.


%% 3. sweep run nmf 
% Run NMF using different k (i.e. number of components) (e.g. 1-10, which takes 20s on my laptop). Plot the reconstruction error as a function of 
% the number of components. 
% Make sure to store the resulting Ws and Hs. They will be used later.
% You may use the matlab function nnmf from the statistics and machine
% learning toolbox
% [W,H,D] = nnmf(data.data_normed,k);





%% 4. visualize W, H, data and the reconstructed data
% Pick a k and visualize W, H. Reconstruct the data from the fitted W and
% H, compare the reconstructed data with the original data.
% You may use the function provided: 
% plot_wh(W,H) and plot_original_reconstructed(X,W,H)
% We recommend using k=4.


%% [IMPORTANT] 5. sort W and H for better interpretations; plot the sorted, W, H, data, and the reconstruction

% You might find the visualization in the previous section to be not so
% illuminating. We shall now sort W and H for better interpretations.

% Hint: for each row (cell) in W, find the component where it participates the
% most, and "assign" the cell to that component. Within each component, sort the cells belonging 
% to the component. 

% p.s. This step will be repeated a few times so we recommend wrapping it
% up inside a function. You may want to store the sorted_index for this and
% the next step.




%% 6. compare different NMF with different ranks; compare the old and new sorting on the new fit
% Repeat step 5, but with a higher number of components (e.g., 5).  
% In addition, we could apply the sorting we obtained from before on this
% fitted result. Compare the two sortings on this fit.



%% 7. compare W with the ground truth cluster
% Now we want to check whether NMF gives us similar grouping of cells as
% given by the clustering done in the original paper. The cluster
% assignment is given in data.clust. 

% hint: You may use B=onehotencode(categorical(data.clust),2) to obtain the
% one-hot encoding (vectors with 0s in all but one element and 1 in that
% element). Then use imagesc to visualize the sorted B. Compare with the
% sorted W.


%% 8. compare with pca
% Fit a PCA and repeat step 5. Compare the result with NMF. You may use the pca
% function in the Statistics and Machine Learning Toolbox.

% hint: notice how coeff and score correspond to W and H (check the sizes)


%% 9. compare the factors across two fits
% Fitting a PCA with 5 PCs, one will get the 4 PCs obtained from fitting a
% PCA with 4 PCs. Let's check whether this holds for NMF. Compare the
% similarity bewteen the Ws obtained from a rank-4 NMF and a rank-5 NMF.
% 
% Hint: you may do a 4 x 5 grid, for each subplot (i,j), do a scatter plot
% of the i, j-th column of the corresponding W from the rank-4 and rank-5 NMF.




%% cutsom helper functions
% 



%% helper functions for plotting

function plot_original_reconstructed(X,W,H)
    X_recon = W * H;
    figure
    tiledlayout(1,2,'Padding','compact')
    clims=[0,20];
    nexttile;
    imagesc(X,clims)
    xlabel('gene')
    ylabel('cell')
    daspect([1 1 1])
    title('original')
    nexttile;
    imagesc(X_recon,clims)
    title('reconstructed')
    xlabel('gene')
    ylabel('cell')
    daspect([1 1 1])
end

function plot_wh(W,H)

    figure('Units','normalized','Position',[0.1,0.8,0.6,0.4]);
    subplot(2,3,[1,4])
    Wclim=[min(min(W)),max(max(W))] * 0.8;
    imagesc(W,Wclim)
    colorbar
    title('W')
    xlabel('factors')
    ylabel('cells')
    subplot(2,3,[2,3])
    Hclim=[min(min(H)),max(max(H))] * 0.8;
    imagesc(H,Hclim)
    colorbar
    title('H')
    xlabel('genes')
    ylabel('factors')
    
end