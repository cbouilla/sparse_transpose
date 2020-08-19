#!/usr/bin/env python3

"""Creates charts (duration, speed up, boxplot...) from a CSV file.

   Copyright (c) 2020 Charles Bouillaguet and Jérôme Bonacchi
"""

from os import mkdir
from sys import argv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fieldnames = ["partitioning", "throughput", "copy_pointers", "buckets"]


def read_csv(filename):
    """Reads `filename` as a CSV file, assuming that the columns of the CSV
    file correspond to `fieldnames`.

    Args:
        filename (str): The CSV file.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing the data in `filename`.
    """
    data = pd.read_csv(filename, header=None, names=fieldnames)
    return data


def plot_mean_duration(data, show=True, save=True, save_path=""):
    """Plots data charts for each algorithm for each matrix.

    Args:
        data (pandas.DataFrame): Data shaped like parallel_means_pivot
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if not show and not save:
        return
    fig, axes = plt.subplots()
    print("Plotting duration...", end=' ', flush=True)
    thread = pd.Series(data.index.values, index=data.index)
    data['partitioning'].div(thread).plot(label="partitioning",
                                          linestyle=':', marker='.', color='b', secondary_y=True)
    data['copy_pointers'].div(thread).plot(
        label="copy_pointers", linestyle=':', marker='.', color='r')
    data['buckets'].div(thread).plot(
        label="buckets", linestyle=':', marker='.', color='g')
    axes.set(xlabel='threads', ylabel="mean duration per thread (s)", xlim=(0, len(data) + 1),
             xticks=np.concatenate(([1], 4*np.arange(1, (len(data) + 1) / 4))),
             title="pre_transpose41")
    axes.grid(linestyle='-.')
    axes.legend()
    if save:
        output = save_path + "_mean_duration.svg"
        fig.savefig(output)
    if show:
        plt.show()
    plt.close(fig)
    print("done.")


def plot_duration(data, show=True, save=True, save_path=""):
    """Plots data charts for each algorithm for each matrix.

    Args:
        data (pandas.DataFrame): Data shaped like parallel_means_pivot
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if not show and not save:
        return
    fig, axes = plt.subplots()
    print("Plotting sum of duration...", end=' ', flush=True)
    data['partitioning'].plot(label="partitioning",
                              linestyle=':', marker='.', color='b', secondary_y=True)
    data['copy_pointers'].plot(
        label="copy_pointers", linestyle=':', marker='.', color='r')
    data['buckets'].plot(label="buckets", linestyle=':', marker='.', color='g')
    axes.set(xlabel='threads', ylabel="duration (s) (sum for each thread)", xlim=(0, len(data) + 1),
             xticks=np.concatenate(([1], 4*np.arange(1, (len(data) + 1) / 4))),
             title="pre_transpose41")
    axes.grid(linestyle='-.')
    axes.legend()
    if save:
        output = save_path + "_duration.svg"
        fig.savefig(output)
    if show:
        plt.show()
    plt.close(fig)
    print("done.")


def plot_throughput(data, show=True, save=True, save_path=""):
    """Plots data charts for each algorithm for each matrix.

    Args:
        data (pandas.DataFrame): Data shaped like parallel_means_pivot
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if not show and not save:
        return
    fig, axes = plt.subplots()
    print("Plotting throughput...", end=' ', flush=True)
    data['throughput'].plot(linestyle=':', marker='.', color='b')
    axes.set(xlabel='threads', ylabel="throughput (GB/s) (16nnz/time_master)", xlim=(0, len(data) + 1),
             xticks=np.concatenate(([1], 4*np.arange(1, (len(data) + 1) / 4))),
             title="pre_transpose41")
    axes.grid(linestyle='-.')
    if save:
        output = save_path + "_throughput.svg"
        fig.savefig(output)
    if show:
        plt.show()
    plt.close(fig)
    print("done.")


def main():
    filename = "csv/bench.csv"
    save_path = "charts/"
    argc = len(argv)
    if argc == 1:
        print("Input file:", filename)
        start_index = filename.find("/") + 1
        save_path = save_path + filename[start_index:-4]
        print("output file:", save_path + "*.svg")
    elif argc == 2:
        filename = argv[1]
        print("Input file:", filename)
        start_index = filename.find("/") + 1
        save_path = save_path + filename[start_index:-4]
        print("output file:", save_path + "*.svg")
    else:
        print(
            "Usage: plot output_filename.")
        return
    data = read_csv(filename)
    data.index = data.index + 1
    plot_mean_duration(data, show=False, save=True, save_path=save_path)
    plot_duration(data, show=False, save=True, save_path=save_path)
    plot_throughput(data, show=False, save=True, save_path=save_path)


if __name__ == "__main__":
    main()
