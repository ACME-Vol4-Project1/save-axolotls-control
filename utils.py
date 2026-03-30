import numpy as np
from matplotlib import pyplot as plt

# TODO: Add easy plotting utilities to produce plots like those in test_model.ipynb
def popplot(sol, title="title", vert=True):
    # plots
    if vert:
        fig, ax = plt.subplots(2, 1)
    else:
        fig, ax = plt.subplots(1, 2)

    y_init = sol.y[:, 0]

    # linear scale
    # graph each compartment
    ax[0].plot(sol.t, sol.y[0] / np.sum(y_init), label='S')
    ax[0].plot(sol.t, sol.y[1] / np.sum(y_init), label="I")
    ax[0].plot(sol.t, sol.y[2] / np.sum(y_init), label="T1")
    ax[0].plot(sol.t, sol.y[3] / np.sum(y_init), label="D")

    # label the graph
    ax[0].set_title("Linear scale")
    ax[0].set_xlabel("time (days)")
    ax[0].set_ylabel("proportion")

    # log log
    # graph each compartment
    ax[1].semilogy(sol.t, sol.y[0] / np.sum(y_init))
    ax[1].semilogy(sol.t, sol.y[1] / np.sum(y_init))
    ax[1].semilogy(sol.t, sol.y[2] / np.sum(y_init))
    ax[1].semilogy(sol.t, sol.y[3] / np.sum(y_init))

    # label the graph
    ax[1].set_title("Log scale")
    ax[1].set_xlabel("time (days)")
    ax[1].set_ylabel("proportion")

    # figure
    fig.legend(bbox_to_anchor=(1.1, 0.6))
    fig.suptitle(title)
    fig.tight_layout()