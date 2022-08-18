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

data.data_normed = data.data ./ std(data.data,[],1);

%% 2. visualize data
% Visualize the normalized data. You might use imagesc.
figure
imagesc(data.data_normed)
colorbar

%% 3. sweep run nmf 
% Run NMF using different k(i.e. number of components) (e.g. 1-10, which takes 20s on my laptop). Plot the reconstruction error as a function of 
% the number of components. 
% Make sure to store the resulting Ws and Hs. They will be used later.
% You may use the matlab function nnmf from the statistics and machine
% learning toolbox
% [W,H,D] = nnmf(data.data_normed,k);

rng(2)
k_to_sweep=1:10; % about 20s
error_l = zeros(length(k_to_sweep),1);
H_l = {};
W_l = {};
for kk=1:length(k_to_sweep)
    k = k_to_sweep(kk);
    [W,H,D] = nnmf(data.data_normed,k);
    H_l{kk} = H; 
    W_l{kk} = W;
    error_l(kk)=D;
end
figure;plot(k_to_sweep,error_l)




%% 4. visualize W, H, data and the reconstructed data
% Pick a k and visualize W, H. Reconstruct the data from the fitted W and
% H, compare the reconstructed data with the original data.
% You may use the function provided: 
% plot_wh(W,H) and plot_original_reconstructed(X,W,H)
% We recommend using k=4.

k=4;
W = W_l{k};
H = H_l{k};
plot_wh(W,H)
plot_original_reconstructed(data.data_normed,W,H)

%% 5. sort W and H for better interpretations; plot the sorted, W, H, data, and the reconstruction

% You might find the visualization in the previous section to be not so
% illuminating. We shall now sort W and H for better interpretations.

% Hint: for each row (cell) in W, find the component where it participates the
% most, and "assign" the cell to that component. Within each component, sort the cells belonging 
% to the component. 

% p.s. This step will be repeated a few times so we recommend wrapping it
% up inside a function. You may want to store the sorted_index for this and
% the next step.


W_ind_sorted_l = {};
H_ind_sorted_l = {};
k=4;
W = W_l{k};
H = H_l{k};

[W_ind_sorted_l{k},H_ind_sorted_l{k}]=sort_and_plot_wh_x(data.data_normed,W,H);


%% 6. compare different NMF with different ranks; compare the old and new sorting on the new fit
% Repeat step 5, but with a higher number of components (e.g., 5).  
% In addition, we could apply the sorting we obtained from before on this
% fitted result. Compare the two sortings on this fit.

k=5;
W = W_l{k};
H = H_l{k};
plot_wh(W(W_ind_sorted_l{4},:),H(:,H_ind_sorted_l{4}))
sgtitle("using old sorting")
[W_ind_sorted_l{k},H_ind_sorted_l{k}]=sort_and_plot_wh_x(data.data_normed,W,H);


%% 7. compare W with the ground truth cluster
% Now we want to check whether NMF gives us similar grouping of cells as
% given by the clustering done in the original paper. The cluster
% assignment is given in data.clust. 
% Do the comparison for both n_components=4 and 5.

% hint: You may use B=onehotencode(categorical(data.clust),2) to obtain the
% one-hot encoding (vectors with 0s in all but one element and 1 in that
% element). Then use imagesc to visualize the sorted B. Compare with the
% sorted W.
B = onehotencode(categorical(data.clust),2);

k_to_compare = [4,5];
n_k_to_compare = length(k_to_compare);
figure

for kk=1:n_k_to_compare
    k = k_to_compare(kk);
%     ax=nexttile(tlo);
    
    subplot(n_k_to_compare,2,(kk-1)*2+1)
    imagesc(B(W_ind_sorted_l{k},:))
    xticks(1:size(B,2))
    xtickformat('%d')
    xlabel('cell cluster')
    ylabel('cell')
    title(['cluster, sorted using k=',num2str(k)])
    
    subplot(n_k_to_compare,2,(kk-1)*2+2)
    imagesc(W_l{k}(W_ind_sorted_l{k},:))
    xlabel('factors')
    ylabel('cell')
    title(['W, k=',num2str(k)])

end

%% 8. compare with pca
% Fit a PCA and repeat step 5. Compare the result with NMF. You may use the pca
% function in the Statistics and Machine Learning Toolbox.

% hint: notice how coeff and score correspond to W and H (check the sizes)

k=4;
% [coeff,score,latent] = pca(data.data_normed,'NumComponents',k);


[score_ind_sorted,coeff_ind_sorted]=sort_and_plot_wh_x(data.data_normed,score,coeff');


%% 9. compare the factors across two fits
% Fitting a PCA with 5 PCs, one will get the 4 PCs obtained from fitting a
% PCA with 4 PCs. Let's check whether this holds for NMF. Compare the
% similarity bewteen the Ws obtained from a rank-4 NMF and a rank-5 NMF.
% 
% Hint: you may do a 4 x 5 grid, for each subplot (i,j), do a scatter plot
% of the i, j-th column of the corresponding W from the rank-4 and rank-5 NMF.

k_to_compare = [4,5];

figure
tiledlayout(k_to_compare(1),k_to_compare(2))
for ii=1:k_to_compare(1)
    for jj=1:k_to_compare(2)  
        nexttile
        scatter(W_l{k_to_compare(1)}(:,ii),W_l{k_to_compare(2)}(:,jj))
        xlabel(sprintf('W component %d from k=%d',ii,k_to_compare(1)))
        ylabel(sprintf('W component %d from k=%d',jj,k_to_compare(2)))
    end             
end



%% cutsom helper functions
%
function [W_ind_sorted,H_ind_sorted]=sort_and_plot_wh_x(X,W,H)
    [W_sorted,W_ind_sorted]=sort_factors(W,2);
    [H_sorted,H_ind_sorted]=sort_factors(H,1);
    
    plot_wh(W_sorted,H_sorted)
    X_sorted = X(W_ind_sorted,H_ind_sorted);
    plot_original_reconstructed(X_sorted ,W_sorted,H_sorted)

end

function [X_sorted,ind_sorted]=sort_factors(X,dim)
% dim: the dimension corresponding to the components, e.g. for W: ncells x ncomponents, dim=2

    [M,I] = max(X,[],dim); % get the max value, as well as which component achieves the max, for each cell (row in W) / gene (col in H)
    ind_sorted=[]; 
    unique_l = unique(I);

    if dim==1 % tranpose to make the code below consistent for both W and H
        M = M';
        I = I';
    end
    
    for i=1:length(unique_l) 
        mask = I == i; % get the boolean mask for which cells/genes belong to component i
        mask_inds = find(mask); % get the indices for cells/genes that belong to component i
        [~,ind_sorted_within]=sort(M(mask)); % sort the max value of the cells/genes within component i
        ind_sorted =[ind_sorted; mask_inds(ind_sorted_within)]; % append the index to the result

    end
    if dim==1
        X_sorted=X(:,ind_sorted);
    else
        X_sorted =X(ind_sorted,:);
    end


end



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