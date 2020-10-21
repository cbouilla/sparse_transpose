#!/usr/bin/env python3

"""Creates charts (duration, speed up, boxplot...) from a CSV file.

   Copyright (c) 2020 Charles Bouillaguet and Jérôme Bonacchi
"""

from os import mkdir
from sys import argv

import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Liberation Serif'],
    'font.size': 9
})
import numpy as np
import pandas as pd
from matplotlib._color_data import CSS4_COLORS


fieldnames = ["program", "cflags", "cxxflags", "threads", "matrix", "compress",
              "transpose", "compress_tr", "transpose_tr", "step", "total_step"]


def extract_matrix_name(path):
    """Extracts the name of the matrix from `path`.

    Args:
        path (str): The path to the matrix file.

    Returns:
        str: The name of the matrix.
    """
    begin = path.rfind('/')
    matrix_name = path[begin + 1:-4]
    return matrix_name


def extract_program_name(name):
    """Extracts the name of the program.

    Args:
        name (str): The raw name of the program.

    Returns:
        str: The correct name of the program.
    """
    program_name = name.replace('::', ' ')
    return program_name


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


def is_serial(name):
    """Returns whether or not algorithm `name` corresponds to a serial
    algorithm.

    Args:
        name (str): The name of the algorithm.

    Returns:
        bool: Whether or not algorithm `name` corresponds to a serial
        algorithm.
    """
    is_gustavson = name.find("Gustavson") != -1
    is_std_sort = name .find("STL") != -1
    is_mkl_serial = (name.find("MKL") != -1) and (name.find("serial") != -1)
    isSerial = is_gustavson or is_std_sort or is_mkl_serial
    return isSerial


def means(data):
    """Computes the means of the several repetitions of the executions of each
    algorithm in `data`.

    Args:
        data (pandas.DataFrame): Data shaped like the output of `read_csv`.

    Returns:
        pandas.DataFrame: The several repetitions of the executions of each
    algorithm are replaced by the mean of these repetitions.
    """
    groups = ["program", "cflags", "cxxflags", "matrix", "threads"]
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
    sequential_indices = (data["program"]).apply(is_serial)
    sequential = data.loc[sequential_indices]
    parallel = data.loc[sequential_indices == False]
    return sequential, parallel


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
    total = dict()
    threads = data.iloc[0]["transpose"].index
    length = len(threads)
    for index in data.index:
        name, cflags, cxxflags, matrix = index
        name = extract_program_name(name)
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
        mean =  data.loc[index, "transpose mean"]
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
    for (k, d) in total.items():
        name, _, _ = k
        fig, axes = plt.subplots()
        transpose, transpose_tr = d['transpose'], d['transpose_tr']
        print("Plotting overall duration of {}...".format(name), end=' ',
              flush=True)
        transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
        transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
                          color='r')
        mean =  transpose.add(transpose_tr).div(2)
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
        data (pandas.DataFrame): Data shaped like parallel_means_pivot
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
        name = extract_program_name(name)
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
        # fig, axes = plt.subplots()
        # print("Plotting speed up of {} on {}...".format(name, matrix), end=' ',
        #       flush=True)
        sequential = ref.loc[ref['matrix']
                             == matrix, 'transpose'].values[0]
        total[key]['sequential'] = total[key]['sequential'] + sequential
        total[key]['overall_transpose'] = total[key]['overall_transpose'].add(
            data.loc[index, "transpose"])
        transpose = data.loc[index, "transpose"].rdiv(sequential)
        total[key]['mean_transpose'] = total[key]['mean_transpose'].add(
            transpose)
#        transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
        sequential_tr = ref.loc[ref['matrix'] ==
                                matrix, 'transpose_tr'].values[0]
        total[key]['sequential_tr'] = total[key]['sequential_tr'] + sequential_tr
        total[key]['overall_transpose_tr'] = total[key]['overall_transpose_tr'].add(
            data.loc[index, "transpose_tr"])
        transpose_tr = data.loc[index, "transpose_tr"].rdiv(sequential_tr)
        total[key]['mean_transpose_tr'] = total[key]['mean_transpose_tr'].add(
            transpose_tr)
        # transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
        #                   color='r')
        # mean = transpose.add(transpose_tr).div(2)
        # mean.plot(label="mean", linestyle=':', marker='.', color='g')
        # axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
        #          xticks=np.concatenate(
        #                 ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
        #          ylabel='speed up', title=name + " on " + matrix)
        # axes.grid(linestyle='-.')
        # axes.legend()
        # if save:
        #     output = save_path + \
        #         name.replace(" ", "_") + "_" + matrix + "_speed_up.svg"
        #     fig.savefig(output)
        # if show:
        #     plt.show()
        # plt.close(fig)
        # print("done.")
    # Mean
    # for (k, d) in total.items():
    #     name, _, _ = k
    #     length = len(total[k]['matrices'])
    #     fig, axes = plt.subplots()
    #     transpose = d['mean_transpose'].div(length)
    #     transpose_tr = d['mean_transpose_tr'].div(length)
    #     print("Plotting mean speed up of {}...".format(name), end=' ',
    #           flush=True)
    #     transpose.plot(label="transpose", linestyle=':', marker='.', color="b")
    #     transpose_tr.plot(label="transpose_tr", linestyle=':', marker='.',
    #                       color='r')
    #     mean = transpose.add(transpose_tr).div(2)
    #     mean.plot(label="mean", linestyle=':', marker='.', color='g')
    #     axes.set(xlabel='threads', xlim=(threads[0] - 1, threads[-1] + 1),
    #              xticks=np.concatenate(
    #                     ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
    #              ylabel='speed up', title=name + " mean speed up")
    #     axes.grid(linestyle='-.')
    #     axes.legend()
    #     if save:
    #         output = save_path + \
    #             name.replace(" ", "_") + "_mean_speed_up.svg"
    #         fig.savefig(output)
    #     if show:
    #         plt.show()
    #     plt.close(fig)
    #     print("done.")
    # Overall
    for (k, d) in total.items():
        name, _, _ = k
        fig, axes = plt.subplots(figsize=(5.5,3.7))
        transpose = d['overall_transpose'].rdiv(d['sequential'])
        transpose_tr = d['overall_transpose_tr'].rdiv(
            d['sequential_tr'])
        print("Plotting overall speed up of {}...".format(name), end=' ',
              flush=True)
        transpose.plot(label="$A^T$", linestyle='', marker='.', markersize=3, color="b")
        transpose_tr.plot(label="$(A^T)^T$", linestyle='', marker='.', markersize=3,
                          color='r')
        mean = transpose.add(transpose_tr).div(2)
        mean.plot(label="moyenne", linestyle='', marker='.', markersize=3, color='g')
        axes.set(xlabel="Fils d'exécution", xlim=(threads[0] - 1, threads[-1] + 1),
                 xticks=np.concatenate(
                        ([1], 4*np.arange(threads[0], (threads[-1] + 1) / 4))),
                 ylabel='Accélération')
        axes.grid(linestyle='-.')
        leg = axes.legend(ncol=3)
#        leg.get_frame().set_alpha(0.6)
        if save:
            output = save_path + \
                name.replace(" ", "_") + "_overall_speed_up.eps"
            fig.savefig(output, dpi=1200)
        if show:
            plt.show()
        plt.close(fig)
        print("done.")


def find_minima(data, parallel):
    """Finds minima (according to the number of threads) of each execution of
    each algorithm in `data` and returns a boolean vector where the executions
    of the minima are in `parallel`.

    Args:
        data (pandas.DataFrame): Data shaped like parallel_means_pivot
        parallel: Data shaped like the output of `read_csv`.

    Returns:
        pandas.Series: A vector with the same indices as in `parallel. and
        boolean values corresponding to whether or not this index is the best
        execution (according to the number of threads).
    """
    parallel_idx = pd.Series([False] * len(parallel),
                             index=parallel.index, dtype="boolean",
                             name="program")
    for data_idx in data.index:
        program, cflags, cxxflags, matrix = data_idx
        thread_min = data.loc[data_idx, "transpose mean"].idxmin()
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
        name = columns_names[0].replace(" ", "\n", 1)
        name = extract_program_name(name)
        labels.append(name)
    x = range(1, len(labels) + 1)
    total = [0] * len(x)
    total_tr = [0] * len(x)
    for i in range(0, len(data.index), N_repeat):
        matrix, _ = data.iloc[i].name
        if matrix.find("pre_transpose6") < 0:
            continue
        # dropna is needed because other threads appear as NaN
        current_data = data[i:i + N_repeat].dropna(axis=1)
        fig, axes = plt.subplots(figsize=(7, 4))
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
        axes.plot(x, means, linestyle="", color='g', marker='.', markersize=3, zorder=10,
                  label='moyenne des medianes')
        axes.set(xticks=x, ylabel='Durée (s)')
        bplot['medians'][0].set_label("$A^T$")
        bplot_tr['medians'][0].set_label("$(A^T)^T$")
        labels = list()
        columns = current_data['transpose'].columns
        for columns_names in columns:
            name = columns_names[0].replace(" ", "\n", 1)
            name = extract_program_name(name)
#            if not is_serial(name):
#                name = name + " (" + str(columns_names[3]) + ")"
            labels.append(name)
        plt.setp(axes, xticklabels=labels)
        plt.xticks(horizontalalignment='center')
        axes.grid(False)
        axes.grid(linestyle='-.', axis='y')
        axes.legend(ncol=3)
        if save:
            output = save_path + matrix + "_boxplot.eps"
            fig.savefig(output, dpi=1200)
        if show:
            plt.show()
        plt.close(fig)
        print("done.")
    # Overall
    # fig, axes = plt.subplots(figsize=(7,4))
    # print("Plotting overall boxes...", end=' ', flush=True)
    # axes.plot(x, total, label="$A^T$", linestyle="", color="b", marker='+', markersize=3)
    # axes.plot(x, total_tr, label="$(A^T)^T$", linestyle="", color="r",
    #           marker='+', markersize=3)
    # means = np.add(total, total_tr) / 2
    # axes.plot(x, means, label="moyenne des médianes", linestyle="", color="g", marker='.', markersize=3)
    # axes.set(xticks=x, ylabel='Durée (s)')
    # plt.setp(axes, xticklabels=labels)
    # plt.xticks(horizontalalignment='center')
    # axes.grid(False)
    # axes.grid(linestyle='-.', axis='y')
    # axes.legend(ncol=3)
    # if save:
    #     output = save_path + "overall_suitesparse_boxplot.eps"
    #     fig.savefig(output, dpi=1200)
    # if show:
    #     plt.show()
    # plt.close(fig)
    # print("done.")
    # print(total)
    # print(total_tr)
    # print(means)


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
        return
    data = read_csv(filename)
    # sequential, parallel = split(data)
    # parallel_means = means(parallel.iloc[:, :-2])
    # transpose_means = parallel_means[[
    #     'transpose', 'transpose_tr']].apply(np.mean, axis=1)
    # parallel_means.insert(loc=len(parallel_means.columns), column="transpose mean",
    #                       value=transpose_means)
    # ref = read_csv("csv/classical_O2.csv")
    # ref_means = means(ref)
    # parallel_means_pivot = parallel_means.pivot_table(
    #     index=['program', 'cflags', 'cxxflags', 'matrix'], columns="threads",
    #     values=['compress', 'transpose', 'compress_tr',
    #             'transpose_tr', 'transpose mean'],
    #     observed=True)
    # # plot_duration(parallel_means_pivot, show=False,
    # #            save=True, save_path=save_path)
    # plot_speed_up(parallel_means_pivot, ref_means,
    #               show=False, save=True, save_path=save_path)
    # best = find_minima(parallel_means_pivot, parallel)
    # faster_parallel = parallel.loc[best]
    # new_data = pd.concat([sequential, faster_parallel])
    new_data = data.query('threads == 1')
    new_data_pivot = new_data.pivot_table(
        index=['matrix', 'step'],
        columns=['program', 'cflags', 'cxxflags', 'threads', 'total_step'],
        values=['compress', 'transpose', 'compress_tr', 'transpose_tr'],
        observed=True)
    plot_box(new_data_pivot, show=False, save=True, save_path=save_path)


if __name__ == "__main__":
    main()
