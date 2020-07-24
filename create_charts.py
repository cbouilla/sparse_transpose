#!/usr/bin/env python3

"""Creates charts (duration, speed up, boxplot...) from a CSV file.

   Copyright (c) Charles Bouillaguet and Jérôme Bonacchi
"""

from os import mkdir
from sys import argv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib._color_data import CSS4_COLORS


fieldnames = ["program", "cflags", "cxxflags", "threads", "matrix", "compress",
              "transpose", "compress_tr", "transpose_tr", "step", "total_step"]


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
    if not show and not save:
        return
    total = dict()
    threads = data.iloc[0]["transpose"].index
    length = len(threads)
    for index in data.index:
        name, cflags, cxxflags, matrix = index
        key = (name, cflags, cxxflags)
        if key not in total:
            total[key] = dict(transpose=pd.Series([0] * length, index=threads),
                              transpose_tr=pd.Series([0] * length, index=threads))
        fig, axes = plt.subplots()
        print("Plotting duration of {} on {}...".format(name, matrix), end=' ',
              flush=True)
        transpose = data.loc[index, "transpose"]
        total[key]['transpose'] = total[key]['transpose'].add(transpose)
        transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
        transpose_tr = data.loc[index, "transpose_tr"]
        total[key]['transpose_tr'] = total[key]['transpose_tr'].add(
            transpose_tr)
        transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
                          color='r')
        mean = transpose.add(transpose_tr).div(2)
        mean.plot(label="mean", linestyle=':', marker='.', color='g')
        axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
                 xticks=np.concatenate(
                        ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
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
    # Overall
    for d in total.values():
        fig, axes = plt.subplots()
        transpose, transpose_tr = d['transpose'], d['transpose_tr']
        print("Plotting overall duration of {}...".format(name), end=' ',
              flush=True)
        transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
        transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
                          color='r')
        mean = transpose.add(transpose_tr).div(2)
        mean.plot(label="mean", linestyle=':', marker='.', color='g')
        axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
                 xticks=np.concatenate(
                        ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
                 ylabel='duration (s)', title=name + " overall")
        axes.grid(linestyle='-.')
        axes.legend()
        if save:
            output = save_path + \
                name.replace(" ", "_") + "_overall_duration.svg"
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
    if not show and not save:
        return
    total = dict()
    threads = data.iloc[0]["transpose"].index
    length = len(threads)
    for index in data.index:
        name, cflags, cxxflags, matrix = index
        key = (name, cflags, cxxflags)
        if key not in total:
            total[key] = dict(sequential=0, sequential_tr=0, matrices=set(),
                              mean_transpose=pd.Series([0] * length,
                                                       index=threads),
                              mean_transpose_tr=pd.Series([0] * length,
                                                          index=threads),
                              overall_transpose=pd.Series([0] * length,
                                                          index=threads),
                              overall_transpose_tr=pd.Series([0] * length,
                                                             index=threads))
        try:
            total[key]['matrices'].add(matrix)
        except:
            pass
        fig, axes = plt.subplots()
        print("Plotting speed up of {} on {}...".format(name, matrix), end=' ',
              flush=True)
        sequential = ref.loc[ref['matrix']
                             == matrix, 'transpose'].values[0]
        total[key]['sequential'] = total[key]['sequential'] + sequential
        total[key]['overall_transpose'] = total[key]['overall_transpose'].add(
            data.loc[index, "transpose"])
        transpose = data.loc[index, "transpose"].rdiv(sequential)
        total[key]['mean_transpose'] = total[key]['mean_transpose'].add(
            transpose)
        transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
        sequential_tr = ref.loc[ref['matrix'] ==
                                matrix, 'transpose_tr'].values[0]
        total[key]['sequential_tr'] = total[key]['sequential_tr'] + sequential_tr
        total[key]['overall_transpose_tr'] = total[key]['overall_transpose_tr'].add(
            data.loc[index, "transpose_tr"])
        transpose_tr = data.loc[index, "transpose_tr"].rdiv(sequential_tr)
        total[key]['mean_transpose_tr'] = total[key]['mean_transpose_tr'].add(
            transpose_tr)
        transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
                          color='r')
        mean = transpose.add(transpose_tr).div(2)
        mean.plot(label="mean", linestyle=':', marker='.', color='g')
        axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
                 xticks=np.concatenate(
                        ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
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
    # Mean
    for k, d in total.items():
        name, _, _ = k
        length = len(total[k]['matrices'])
        fig, axes = plt.subplots()
        transpose = d['mean_transpose'].div(length)
        transpose_tr = d['mean_transpose_tr'].div(length)
        print("Plotting mean speed up of {}...".format(name), end=' ',
              flush=True)
        transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
        transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
                          color='r')
        mean = transpose.add(transpose_tr).div(2)
        mean.plot(label="mean", linestyle=':', marker='.', color='g')
        axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
                 xticks=np.concatenate(
                        ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
                 ylabel='speed up', title=name + " mean speed up")
        axes.grid(linestyle='-.')
        axes.legend()
        if save:
            output = save_path + \
                name.replace(" ", "_") + "_mean_speed_up.svg"
            fig.savefig(output)
        if show:
            plt.show()
        plt.close(fig)
        print("done.")
    # Overall
    for k, d in total.items():
        name, _, _ = k
        fig, axes = plt.subplots()
        transpose = d['overall_transpose'].rdiv(d['sequential'])
        transpose_tr = d['overall_transpose_tr'].rdiv(
            d['sequential_tr'])
        print("Plotting overall speed up of {}...".format(name), end=' ',
              flush=True)
        transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
        transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
                          color='r')
        mean = transpose.add(transpose_tr).div(2)
        mean.plot(label="mean", linestyle=':', marker='.', color='g')
        axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
                 xticks=np.concatenate(
                        ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
                 ylabel='speed up', title=name + " overall speed up")
        axes.grid(linestyle='-.')
        axes.legend()
        if save:
            output = save_path + \
                name.replace(" ", "_") + "_overall_speed_up.svg"
            fig.savefig(output)
        if show:
            plt.show()
        plt.close(fig)
        print("done.")
    print(total)


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
                             index=parallel.index, dtype="boolean",
                             name="program")
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
    if not show and not save:
        return
    N_repeat = data.columns[0][-1]
    labels = list()
    columns = data[:1]['transpose'].dropna(axis=1).columns
    for columns_names in columns:
        labels.append(columns_names[0].replace(" ", "\n", 1))
    x = range(1, len(labels) + 1)
    total = [0] * len(x)
    total_tr = [0] * len(x)
    # for columns_names in current_data['transpose_tr'].columns:
    # labels.append(columns_names[0].replace(" ", "\n", 1))
    # labels = labels * 2
    for i in range(0, len(data.index), N_repeat):
        matrix, _ = data.iloc[i].name
        # dropna is needed because other threads appear as NaN
        current_data = data[i:i + N_repeat].dropna(axis=1)
        fig, axes = plt.subplots()
        print("Plotting boxes on {}...".format(matrix), end=' ', flush=True)
        color = {'boxes': CSS4_COLORS['dimgrey'], 
                 'whiskers': CSS4_COLORS['dimgrey'],
                 'medians': 'blue', 'caps': CSS4_COLORS['royalblue']}
        bplot = current_data["transpose"].boxplot(
            color=color, sym="b+", return_type="dict")
        color = {'boxes': CSS4_COLORS['dimgrey'],
                 'whiskers': CSS4_COLORS['dimgrey'],
                 'medians': 'red', 'caps': CSS4_COLORS['firebrick']}
        bplot_tr = current_data["transpose_tr"].boxplot(
            color=color, sym="r+", return_type="dict")
        means = [0] * len(x)
        for i, (median, median_tr) in enumerate(zip(bplot['medians'], 
                                                    bplot_tr['medians'])):
            y = median.get_ydata()[0]
            y_tr = median_tr.get_ydata()[0]
            means[i] = (y + y_tr) / 2
            total[i] = total[i] + y
            total_tr[i] = total_tr[i] + y_tr
        axes.plot(x, means, linestyle="", color='g', marker='.', zorder=10,
                  label='mean of medians')
        axes.set(xticks=x, xlabel='name', ylabel='duration (s)', title=matrix)
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
    print("Plotting overall boxes...", end=' ', flush=True)
    axes.plot(x, total, label="transpose", linestyle="", color="b", marker='+')
    axes.plot(x, total_tr, label="transpose_tr", linestyle="", color="r",
              marker='+')
    means = np.add(total, total_tr) / 2
    axes.plot(x, means, label="mean", linestyle="", color="g", marker='+')
    axes.set(xticks=x, xlabel='name', ylabel='duration (s)',
             title="Overall sum of medians")
    plt.setp(axes, xticklabels=labels)
    plt.xticks(horizontalalignment='center')
    axes.grid(False)
    axes.grid(linestyle='-.', axis='y')
    axes.legend()
    if save:
        output = save_path + "overall_boxplot.svg"
        fig.savefig(output)
    if show:
        plt.show()
    plt.close(fig)
    print("done.")


def main():
    filename = "csv/benchmarks.csv"
    save_path = "charts/"
    argc = len(argv)
    if argc == 1:
        print("Input file:", filename)
        print("Output directory:", save_path)
    elif argc == 2:
        filename = argv[1]
        print("Input file:", filename)
        start_index = filename.find("/") + 1
        save_path = save_path + filename[start_index:-4] + "/"
        created = False
        while not created:
            try:
                mkdir(save_path)
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
                  save=False, save_path=save_path)
    plot_speed_up(parallel_means_pivot, ref_means,
                  show=False, save=False, save_path=save_path)
    best = find_minima(parallel_means_pivot, parallel)
    faster_parallel = parallel.loc[best]
    new_data = pd.concat([sequential, faster_parallel])
    new_data_pivot = new_data.pivot_table(
        index=['matrix', 'step'],
        columns=['program', 'cflags', 'cxxflags', 'threads', 'total_step'],
        values=['compress', 'transpose', 'compress_tr', 'transpose_tr'],
        observed=True)
    plot_box(new_data_pivot, show=True, save=False, save_path=save_path)


if __name__ == "__main__":
    main()
