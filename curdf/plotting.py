import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_rdf(bins, gr, path=None, show=False, xlabel="r (Å)", ylabel="g(r)", title=None):
    """
    Plot g(r) versus r and optionally save the figure.

    Parameters
    ----------
    bins
        Bin centers array.
    gr
        g(r) values matching ``bins``.
    path
        Optional output path to save the figure.
    show
        Whether to display the plot interactively.
    xlabel
        X-axis label text.
    ylabel
        Y-axis label text.
    title
        Optional plot title.

    Returns
    -------
    matplotlib.figure.Figure
        Matplotlib figure containing the plot.
    """
    fig, ax = plt.subplots()
    ax.plot(bins, gr)
    ax.axhline(1.0, color="k", alpha=0.5, zorder=0, linestyle="-")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if path:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return fig
