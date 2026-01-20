class SparseMatrix:
    import pandas as pd
    def __init__(self, df:pd.DataFrame, direction:str):
        self.df = df
        self.direction = direction
        self.tfs = []
    
    def add_TF(self, tf_motif):
        """
        Annotate the SparseMatrix with a given transcription factor motif location.
        Args:
            tf_motif (tuple): A tuple in the format of (chromosome, start, end, TF_name) specifying the transcription factor motif region.
        Returns:
            None: The method modifies the SparseMatrix in place by adding a new attribute 'tfs'.
        """
        self.tfs.append(tf_motif)
        return
    
    def cluster_by_TF(self, tf_motif):
        """
        Cluster the SparseMatrix by a given transcription factor motif removing rows that do not overlap the motif. Then performs clustering on the remaining rows.
        Args:
            tf_motif (tuple): A tuple in the format of (chromosome, start, end, TF_name) specifying the transcription factor motif region.
        Returns:
            SparseMatrix: A new SparseMatrix instance containing only the clustered data overlapping the TF motif.
        """
        columns = [col for col in self.df.columns if tf_motif[1] <= col <= tf_motif[2]]

    def plot_region(self, TSS: bool=False, save_as: str=None):
        """
        Plot the SparseMatrix data for visualization.
        Args:
            TSS (bool): Whether to plot the TSS region.
            save_as (str): The filename to save the plot as.
        Returns:
            None: The method generates a plot of the SparseMatrix data.
        """
        import matplotlib.pyplot as plt

        plt.figure(figsize=(10, 6))
        
