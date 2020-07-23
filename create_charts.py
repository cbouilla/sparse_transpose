#!/usr/bin/env python3

"""Creates charts (duration, speed up, boxplot...) from a CSV file.
"""

import matplotlib._color_data as mcd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

# plt.rcParams.update({'font.size': 12})

fieldnames = ["program", "cflags", "cxxflags", "threads", "matrix", "compress",
              "transpose", "compress_tr", "transpose_tr", "step", "total_step"]

# matrices = ["language", "ASIC_680k", "circuit5M", "flickr", "memchip", "rajat21", "sme3Dc", "stomach", "transient", "webbase-1M", "wiki-Talk", "cage14", "eu-2005",
#             "FullChip", "mac_econ_fwd500", "para-4", "rajat29", "Stanford_Berkeley", "torso1", "venkat01", "web-Google"] + ["pre_transpose"+str(i) for i in range(1, 59)]

# names = ["Gustavson", "std::sort", "tbb::sort",
#          "MKL", "ScanTrans", "MergeTrans"]


def extract_matrix_name(path):
    """Extracts the matrix name from `path`.

    Args:
        path (str): Tha path to the matrix file.

    Returns:
        str: The name of the matrix.
    """
    begin = path.rfind('/')
    matrix_name = path[begin + 1:-4]
    return matrix_name


def read_csv(filename):
    """Reads `filename` as a CSV file, assuming that the columns of the CSV
    file correspond to `fieldnames`.

    Args:
        filename (str): The CSV file.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing the data in `filename`.
    """
    data = pd.read_csv(filename, header=None, names=fieldnames,
                       dtype={'threads': "Int8", 'step': "Int8",
                              'total_step': "Int8"},
                       converters={'matrix': extract_matrix_name})
    data["program"] = data["program"].astype('category')
    data['matrix'] = data['matrix'].astype('category')
    return data


def is_sequential(name):
    """Returns whether or not algorithm `name` corresponds to a sequential
    algorithm.

    Args:
        name (str): The name of the algorithm.

    Returns:
        bool: Whether or not algorithm `name` corresponds to a sequential
        algorithm.
    """
    is_gustavson = name.find("Gustavson") != -1
    is_std_sort = name .find("std::sort") != -1
    is_mkl_sequential = (name.find("MKL") != -
                         1) and (name.find("sequential") != -1)
    is_sequential = is_gustavson or is_std_sort or is_mkl_sequential
    return is_sequential


def means(data):
    """Computes the means of the several repetitions of the executions of each
    algorithm in `data`.

    Args:
        data (pandas.DataFrame): Data shaped like the output of `read_csv`.

    Returns:
        pandas.DataFrame: The several repetitions of the executions of each
    algorithm are replaced by the mean of these repetitions.
    """
    groups = ["program", "cflags", "cxxflags", "threads", "matrix"]
    # dropna is needed with categorical data
    means = (data.groupby(groups, as_index=False)).mean().dropna()
    return means


def split(data):
    """Splits `data` between those from sequential algorithms and those from
    parallel algorithms.

    Args:
        data (pandas.DataFrame): Data shaped like the output of `csv_read`.

    Returns:
        (pandas.DataFrame, pandas.DataFrame): The first DataFrame contains the
        elements of `data` where algorithms are sequential whereas the second
        contains those where algorithms are parallel.
    """
    sequential_indices = (data["program"]).apply(is_sequential)
    sequential = data.loc[sequential_indices]
    parallel = data.loc[sequential_indices == False]
    return sequential, parallel


def plot_duration(data, show=True, save=True, save_path=""):
    """Plots data charts for each algorithm for each matrix.

    Args:
        data (pandas.DataFrame): Data shaped like the output of #TODO
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if show or save:
        for index in data.index:
            name, _, _, matrix = index
            fig, axes = plt.subplots()
            print("Plotting duration of", name, "on", matrix, end=' ', flush=True)
            threads = data.loc[index, "transpose"].index
            transpose = data.loc[index, "transpose"]
            transpose.plot(label="transpose", linestyle=':',
                           marker='.', color="b")
            transpose_tr = data.loc[index, "transpose_tr"]
            transpose_tr.plot(label="transpose_tr", linestyle=':',
                              marker='.', color='r')
            mean = transpose.add(transpose_tr).div(2)
            mean.plot(label="mean", linestyle=':', marker='.',
                      color='g')
            axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
                     xticks=np.concatenate(([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
                     ylabel='duration (s)', title=name + " on " + matrix)
            axes.grid(linestyle='-.')
            axes.legend()
            if save:
                output = save_path + \
                    name.replace(" ", "_") + "_" + matrix + "_duration.svg"
                fig.savefig(output)
            if show:
                plt.show()
            plt.close(fig)
            print("done.")


def plot_speed_up(data, ref, show=True, save=True, save_path=""):
    """Plots speed up charts for each algorithm for each matrix.

    Args:
        data (pandas.DataFrame): Data shaped like the output of #TODO
        ref (pandas.DataFrame): Data shaped like the output of `means`.
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if show or save:
        for index in data.index:
            name, _, _, matrix = index
            fig, axes = plt.subplots()
            print("Plotting speed up of", name, "on", matrix, end=' ', flush=True)
            threads = data.loc[index, "transpose"].index
            sequential = ref.loc[ref['matrix']
                                 == matrix, 'transpose'].values[0]
            transpose = data.loc[index, "transpose"].rdiv(sequential)
            transpose.plot(label="transpose", linestyle=':',
                           marker='.', color="b")
            sequential_tr = ref.loc[ref['matrix'] ==
                                    matrix, 'transpose_tr'].values[0]
            transpose_tr = data.loc[index, "transpose_tr"].rdiv(sequential_tr)
            transpose_tr.plot(label="transpose_tr", linestyle=':',
                              marker='.', color='r')
            mean = transpose.add(transpose_tr).div(2)
            mean.plot(label="mean", linestyle=':',
                      marker='.', color='g')
            axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
                     xticks=np.concatenate(([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
                     ylabel='speed up', title=name + " on " + matrix)
            axes.grid(linestyle='-.')
            axes.legend()
            if save:
                output = save_path + \
                    name.replace(" ", "_") + "_" + matrix + "_speed_up.svg"
                fig.savefig(output)
            if show:
                plt.show()
            plt.close(fig)
            print("done.")


def find_minima(data, parallel):
    """Finds minima (according to the number of threads) of each execution of
    each algorithm in `data` and returns a boolean vector where the executions
    of the minima are in `parallel`.

    Args:
        data (pandas.DataFrame): Data shaped like the output of #TODO
        parallel: Data shaped like the output of `read_csv`.

    Returns:
        pandas.Series: A vector with the same indices as in `parallel. and
        boolean values corresponding to whether or not this index is the best
        execution (according to the number of threads).
    """
    parallel_idx = pd.Series([False]*parallel.shape[0],
                             index=parallel.index, dtype="boolean", name="program")
    for data_idx in data.index:
        program, cflags, cxxflags, matrix = data_idx
        thread_min = data.loc[data_idx, "transpose"].idxmin()
        new_idx = (parallel['program'] == program) & (parallel['cflags'] == cflags) & (parallel['cxxflags']
                                                                                       == cxxflags) & (parallel['threads'] == thread_min) & (parallel['matrix'] == matrix)
        parallel_idx = parallel_idx | new_idx
    return parallel_idx


def plot_box(data, show=True, save=True, save_path=""):
    """Plots box for each algorithm for each matrix.

    Args:
        data (pandas.DataFrame): Data shaped like the output of `read_csv`.
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if show or save:
        N_repeat = data.columns[0][-1]
        labels = list()
        columns = data[:1]['transpose'].dropna(axis=1).columns
        for columns_names in columns:
            labels.append(columns_names[0].replace(" ", "\n", 1))
        x = range(1, len(labels) + 1)
        medians = [0] * len(x)
        medians_tr = [0] * len(x)
        # for columns_names in current_data['transpose_tr'].columns:
            # labels.append(columns_names[0].replace(" ", "\n", 1))
        # labels = labels * 2
        for i in range(0, len(data.index), N_repeat):
            matrix, _ = data.iloc[i].name
            # dropna is needed because other threads appear as NaN
            current_data = data[i:i + N_repeat].dropna(axis=1)
            fig, axes = plt.subplots()
            print("Plotting boxes on", matrix, end=' ', flush=True)
            color = {'boxes': mcd.CSS4_COLORS['dimgrey'], 'whiskers': mcd.CSS4_COLORS['dimgrey'],
                     'medians': 'blue', 'caps': mcd.CSS4_COLORS['royalblue']}
            bplot = current_data["transpose"].boxplot(
                color=color, sym="b+", return_type="dict")
            color = {'boxes': mcd.CSS4_COLORS['dimgrey'], 'whiskers': mcd.CSS4_COLORS['dimgrey'],
                     'medians': 'red', 'caps': mcd.CSS4_COLORS['firebrick']}
            bplot_tr = current_data["transpose_tr"].boxplot(
                color=color, sym="r+", return_type="dict")
            means = [0] * len(x)
            for i, (median, median_tr) in enumerate(zip(bplot['medians'], bplot_tr['medians'])):
                y = median.get_ydata()[0]
                y_tr = median_tr.get_ydata()[0]
                means[i] = (y + y_tr) / 2
                medians[i] = medians[i] + y
                medians_tr[i] = medians_tr[i] + y_tr
            axes.plot(x, means, linestyle="", color='g', marker='.', zorder=10,
                      label='mean of medians')
            axes.set(xticks=x, xlabel='name', ylabel='duration (s)',
                     title=matrix)
            bplot['medians'][0].set_label("transpose")
            bplot_tr['medians'][0].set_label("transpose_tr")
            plt.setp(axes, xticklabels=labels)
            plt.xticks(horizontalalignment='center')
            axes.grid(False)
            axes.grid(linestyle='-.', axis='y')
            axes.legend()
            if save:
                output = save_path + matrix + "_boxplot.svg"
                fig.savefig(output)
            if show:
                plt.show()
            plt.close(fig)
            print("done.")
        # Overall
        fig, axes = plt.subplots()
        print("Overall boxplot", end=' ')
        axes.plot(x, medians, label="transpose", linestyle="", color="b", marker='+')
        axes.plot(x, medians_tr, label="transpose_tr", linestyle="", color="r", marker='+')
        axes.plot(x, np.add(medians, medians_tr) / 2, label="mean", linestyle="", color="g", marker='+')
        axes.set(xticks=x, xlabel='name',
                 ylabel='duration (s)',
                 title="Overall sum of medians")
        plt.setp(axes, xticklabels=labels)
        plt.xticks(horizontalalignment='center')
        axes.grid(False)
        axes.grid(linestyle='-.', axis='y')
        axes.legend()
        if save:
            output = save_path + matrix + "_overall.svg"
            fig.savefig(output)
        if show:
            plt.show()
        plt.close(fig)
        print("done.")


def main():
    filename = "csv/benchmarks.csv"
    save_path = "charts/"
    argc = len(sys.argv)
    if argc == 1:
        print("Input file:", filename)
        print("Output directory:", save_path)
    elif argc == 2:
        filename = sys.argv[1]
        print("Input file:", filename)
        start_index = filename.find("/") + 1
        save_path = save_path + filename[start_index:-4] + "/"
        created = False
        while not created:
            try:
                os.mkdir(save_path)
                created = True
            except FileExistsError:
                if save_path == "charts/":
                    created = True
                    break
                else:
                    print("The directory '{}' already exists.".format(save_path))
                    save_path = input("Choose a new output directory: ")
                    if save_path[-1] != "/":
                        save_path = save_path + "/"
        print("Output directory:", save_path)
    elif argc > 2:
        print(
            "Usage: create_charts [input_filename].\nThe default filename is",
            filename)
        # return
    data = read_csv(filename)
    sequential, parallel = split(data)
    parallel_means = means(parallel.iloc[:, :-2])
    ref = read_csv("csv/classical_O2.csv")
    ref_means = means(ref)
    parallel_means_pivot = parallel_means.pivot_table(
        index=['program', 'cflags', 'cxxflags', 'matrix'], columns="threads",
        values=['compress', 'transpose', 'compress_tr', 'transpose_tr'],
        observed=True)
    plot_duration(parallel_means_pivot, show=False,
                    save=True, save_path=save_path)
    plot_speed_up(parallel_means_pivot, ref_means,
                    show=False, save=True, save_path=save_path)
    best = find_minima(parallel_means_pivot, parallel)
    faster_parallel = parallel.loc[best]
    # parallel_means_pivot.iloc[0].min(level=0)
    # parallel_means[parallel_means == parallel_means_pivot.iloc[0].min(level=0)].iloc[:,-4:].dropna()
    new_data = pd.concat([sequential, faster_parallel])
    # new_data_pivot = new_data.pivot_table(index=['program', 'cflags', 'cxxflags', 'threads', 'step', 'total_step'],
    #                                                   columns=['matrix'],
    #                                                   values=['compress', 'transpose', 'compress_tr', 'transpose_tr'])
    # new_data_pivot = new_data.pivot_table(index=['program', 'cflags', 'cxxflags', 'matrix', 'threads', 'total_step', 'step'],
    #                                                   values=['compress', 'transpose', 'compress_tr', 'transpose_tr'])
    new_data_pivot = new_data.pivot_table(index=['matrix', 'step'],
        columns=['program', 'cflags', 'cxxflags', 'threads', 'total_step'],
        values=['compress', 'transpose', 'compress_tr', 'transpose_tr'],
        observed=True)
    plot_box(new_data_pivot, show=False, save=True, save_path=save_path)


if __name__ == "__main__":
    main()
