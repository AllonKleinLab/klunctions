import numpy as np
import scipy.sparse
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def get_marker_gene_scores(counts, genes, label_list, marker_gene_sets):
    ''' Calculate composite marker gene scores for multiple cell types
    
    Each cell type's scores are calculated as the sum of the 
    z-score-normalized expression values of that cell type's 
    marker genes.
        
    Arguments
    - counts : 2-D scipy.sparse matrix or numpy array, shape (n_cells, n_genes)
        Normalized counts matrix
    - genes : 1-D numpy array (str), shape (n_genes,)
        List of gene names corresponding to the columns of `counts`
    - label_list : list (str), length n_celltypes
        Labels for each cell type
    - marker_gene_sets : list, length n_celltypes
        List of gene sets used to define each cell type.
        Entry j is a list of gene names marking cell type j.
        
    Returns
    scores : 2-D numpy array (float), shape (n_cells, n_celltypes)
        (i,j)th entry is cell i's score for cell type j. 
    '''
    
    n_cells = counts.shape[0]
    n_labels = len(label_list)
    
    # Initialize cell type score matrix - each cell is 
    # assigned a score for each cell type based on its
    # relative expression of cell type marker genes.
    
    scores = np.nan * np.ones((n_cells, n_labels))
    
    for iL, label in enumerate(label_list):
        # Find indices of marker genes
        gene_mask = np.in1d(genes, marker_gene_sets[iL])
        
        # Get expression values of marker genes
        if scipy.sparse.issparse(counts):
            marker_counts = counts[:, gene_mask].toarray()
        else:
            marker_counts = counts[:, gene_mask]
            
        # Calculate cell type score: sum of z-score expression of marker genes
        scores[:, iL] = zscore(marker_counts, ddof=1).sum(1)
    
    return scores

def plot_multiple_scores(scores, embedding, label_list, color_list=None, score_percentile_range=(90, 99.8), ax=None, figure_size=(6,8), point_size=15, show_legend=True, legend_font_size=14):
    ''' Plot cell type scores as different colors on a 2-D embedding
    
    Arguments
    - scores : 2-D numpy array (float), shape (n_cells, n_celltypes)
        (i,j)th entry is cell i's score for cell type j.
    - embedding : 2-D numpy array (float), shape (n_cells, 2)
        Coordinates for plotting cells in 2-D
    - label_list : list (str), length n_celltypes
        Labels for each cell type
    - color_list : list, optional (default: None)
        List of RGB tuples to use for coloring each 
        cell type. If provided, len(color_list)=n_celltypes, 
        with each entry being a list [R,G,B].
        If `None`, colors are n_celltypes evenly distributed
        colors from the jet color map.
    - score_percentile_range : tuple (float), length 2, optional (default = (90, 99.8))
        Range for normalizing each cell type's scores. Score j 
        is normalized by subtracting
        `np.percentile(score[:,j], score_percentile_range[0])`
        and then dividing by 
        `np.percentile(score[:,j], score_percentile_range[1])`, 
        then setting a minimum value of 0 and maximum of 1.
    - ax : matplotlib axes handle, optional
        axes for making drawing the plot

    Returns
    - ax : axes handle
    
    '''
    
    if color_list is None:
        cm = plt.cm.jet
        color_list = []
        for i in range(0, cm.N, int(np.ceil(float(cm.N) / len(label_list)))):
            color_list.append(cm(i)[:3])
    elif len(label_list) != len(color_list):
        print('Error: `label_list` must have the same length as `color_list`')
        return

    # For visualization, we want to emphasize only the highest-scoring
    # cells, so the each cell type's scores are truncated to give
    # non-zero expression to the top 10% of cells. We also reduce
    # sensitivity to outliers by setting at upper bound at the 99.8th
    # percentile.
    scores = scores - np.percentile(scores, score_percentile_range[0], axis=0)
    scores = scores / np.percentile(scores, score_percentile_range[1], axis=0)

    # Set minimum score to 0, maximum to 1
    scores[scores < 0] = 0
    scores[scores > 1] = 1

    # Calculate the colors of each cell
    cell_colors = np.zeros((scores.shape[0], 3))
    for iL, label in enumerate(label_list):
        label_color = np.array(color_list[iL])[None, :]
        color_weighted = np.tile(scores[:, [iL]], (1, 3)) * (1 - np.tile(label_color, (scores.shape[0], 1)))
        cell_colors += color_weighted

    cell_colors[cell_colors > 1] = 1
    cell_colors = 1 - cell_colors

    # The "default" color for cells (those that do not highly
    # express any of the tested cell type programs) are gray. 
    # We set the "amount" of gray as (1 - (max score)).
    gray_custom = 1 - scores.max(1)

    # Adjust cell colors by cell-specific gray level
    cell_colors = cell_colors * (1 - 0.15 * gray_custom[:, None])

    # Cells with high cell type scores will be plotted
    # in front of those with low scores.
    plot_order = np.argsort(scores.max(1))

    # Draw the plot
    if ax is None:
        fig,ax = plt.subplots(figsize=figure_size)

    ax.scatter(embedding[plot_order, 0], embedding[plot_order, 1], c=cell_colors[plot_order, :], s=point_size)
    ax.set_axis_off()

    # Build legend labels
    if show_legend:
        legend_labels = []
        for iL, label in enumerate(label_list):
            legend_labels.append(Line2D([], [], color=color_list[iL], 
                                        marker='o', 
                                        markersize=point_size, 
                                        label=label, 
                                        linewidth=0)
                                )

        # Show legend
        leg = ax.legend(handles = legend_labels, fontsize=legend_font_size)


    return ax
