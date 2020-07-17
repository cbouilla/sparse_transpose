#!/usr/bin/env python3

"""Creates charts (boxplot, speed up) from a CSV file.
    """

import csv
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

plt.rcParams.update({'font.size': 12})

fieldnames = ["name", "cflags", "cxxflags", "threads", "matrix", "compress",
              "transpose", "compress_tr", "transpose_tr", "step", "total_step"]
matrices = [
    "language",
    "ASIC_680k",
    "circuit5M",
    "flickr",
    "memchip",
    "rajat21",
    "sme3Dc",
    "stomach",
    "transient",
    "webbase-1M",
    "wiki-Talk",
    "cage14",
    "eu-2005",
    "FullChip",
    "mac_econ_fwd500",
    "para-4",
    "rajat29",
    "Stanford_Berkeley",
    "torso1",
    "venkat01",
    "web-Google"
] + ["pre_transpose"+str(i) for i in range(1, 59)]
names = ["Gustavson", "std::sort", "tbb::sort",
         "MKL", "ScanTrans", "MergeTrans"]


def csv_read(filename):
    """Reads `filename` as a CSV file, assuming that the columns of the CSV file
    correspond to `fieldnames`.

    Args:
        filename (str): The CSV file.

    Returns:
        [(list of dict, set)]: The first variable is all data in `filename`.
        This is a list where every element corresponds to a dictionary
        containing the values of a row of `filename`. The second variable
        represents the different executions of algorithms. This is a set where
        keys are tuples made of (name, cflags, matrix, threads).
    """
    data = list()
    available = set()
    with open(filename, newline='') as csvfile:
        csvReader = csv.DictReader(csvfile, fieldnames=fieldnames)
        for row in csvReader:
            matrix = row['matrix']
            # Remove file path
            start_index = matrix.rfind('/')
            if start_index != -1:
                matrix = matrix[start_index + 1:]
            # Remove file extension
            matrix = matrix[:-4]
            processed_row = {'name': row['name'],
                             'cflags': row['cflags'],
                             'cxxflags': row['cxxflags'],
                             'threads': int(row['threads']),
                             'matrix': matrix,
                             'compress': float(row['compress']),
                             'transpose': float(row['transpose']),
                             'compress_tr': float(row['compress_tr']),
                             'transpose_tr': float(row['transpose_tr']),
                             'step': int(row['step']),
                             'total_step': int(row['total_step'])
                             }
            data.append(processed_row)
            key = (processed_row['name'], processed_row['cflags'],
                   processed_row['matrix'], processed_row['threads'])
            available.add(key)
    return data, available


def compute_means(data):
    """Computes the means of the several repetitions of the executions of each
    algorithm in `data`.

    Args:
        data (list of dict): Data shaped like the output of `csv_read`.

    Returns:
        list of dict: The several repetitions of the executions of each
    algorithm are replaced by the mean of these repetitions.
    """
    N_repeat = int(data[0]['total_step'])
    means = list()
    for i in range(0, len(data), N_repeat):
        duration = [0, 0, 0, 0]
        for j in range(N_repeat):
            duration[0] += data[i + j]['compress']
            duration[1] += data[i + j]['transpose']
            duration[2] += data[i + j]['compress_tr']
            duration[3] += data[i + j]['transpose_tr']
        for j in range(len(duration)):
            duration[j] = duration[j] / N_repeat
        row = {'name': data[i]['name'], 'cflags': data[i]['cflags'],
               'cxxflags': data[i]['cxxflags'], 'threads': data[i]['threads'],
               'matrix': data[i]['matrix'], 'compress': duration[0],
               'transpose': duration[1], 'compress_tr': duration[2],
               'transpose_tr': duration[3]
               }
        means.append(row)
    return means


def is_sequential(name):
    """Returns whether or not algorithm `name` is sequential.

    Args:
        name (str): The name of the algorithm.

    Returns:
        bool: Whether or not algorithm `name` is sequential.
    """
    is_gustavson = name.find("Gustavson") != -1
    is_std_sort = name .find("std::sort") != -1
    is_mkl_sequential = (name.find("MKL") != -
                         1) and (name.find("sequential") != -1)
    is_sequential = is_gustavson or is_std_sort or is_mkl_sequential
    return is_sequential


def split(data):
    """Splits `data` between those from sequential algorithms and those from
    parallel algorithms.

    Args:
        data (list of dict): Data shaped like the output of `csv_read`.

    Returns:
        (list of dict, list of dict): The first variable contains the elements
        of `data` where algorithms are sequential. The second variable contains
        those where algorithms are parallel.
    """
    sequential = list()
    parallel = list()
    for row in data:
        if is_sequential(row['name']):
            sequential.append(row)
        else:
            parallel.append(row)
    return sequential, parallel


def create_parallel_duration(available):
    """Returns a dictionary where keys are parallel algorithms (ignoring the
    number of threads) in `available` and values are empty lists.

    Args:
        available (set): Data shaped like the output of `csv_read`.

    Returns:
        dict: A dictionary where keys are parallel algorithms in `available`
    and values are empty lists.
    """
    duration = dict()
    # value = {'threads': list(), 'compress': list(), 'transpose': list(),
    #                  'compress_tr': list(), 'transpose_tr': list()}
    value = {'threads': list(), 'transpose': list(), 'transpose_tr': list()}
    for key in available:
        name, cflags, matrix, _ = key
        short_key = (name, cflags, matrix)
        if not is_sequential(name):
            duration[short_key] = value
    return duration


def fill_parallel_duration(data_parallel_means, duration):
    """Fills `duration` with the mean time (from `data_parallel_means`) of
    each execution of each algorithm, puts the executions with different number
    of threads together.

    Args:
        data_parallel_means (list of dict): Data shaped like the output of
        `compute_means`.
        duration (dict): Data shaped like the output of
        `create_parallel_duration.
    """
    a = True  # horrible
    for row in data_parallel_means:
        key = (row['name'], row['cflags'], row['matrix'])
        if a:  # horrible
            a = False
            previous_key = key  # horrible
        # print(key[2], duration[key])
        if key != previous_key:  # horrible
            duration[key] = {'threads': list(), 'transpose': list(),
                             'transpose_tr': list()}
        value = duration[key]
        value['threads'].append(row['threads'])
        value['transpose'].append(row['transpose'])
        value['transpose_tr'].append(row['transpose_tr'])
        # if row['name'] not in {"SpanTrans", "MergeTrans"}:
        #     value['compress'].append(row['compress'])
        #     value['compress_tr'].append(row['compress_tr'])
        previous_key = key  # horrible


def plot_duration(duration, show=True, save=True, save_path=""):
    """Plots duration charts for each algorithm for each matrix.

    Args:
        duration (dict): Data shaped like the output of
        `fill_parallel_duration`
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if show or save:
        for execution in duration:
            name, _, matrix = execution
            threads = duration[execution]['threads']
            time_transpose = duration[execution]['transpose']
            time_transpose_tr = duration[execution]['transpose_tr']
            fig, axes = plt.subplots(figsize=(8, 8))
            # , xlim=(threads[0], threads[-1])
            axes.set(xlabel='threads',
                     xticks=range(threads[0], threads[-1] + 1),
                     ylabel='duration (s)', title=name + " on " + matrix)
            axes.grid(True, linestyle='-.')
            line_transpose, = axes.plot(
                threads, time_transpose, linestyle=(0, (5, 7)), marker='o')
            line_transpose_tr, = axes.plot(threads, time_transpose_tr, 'ro--')
            line_transpose.set_label("transpose")
            line_transpose_tr.set_label("transpose_tr")
            axes.legend()
            if save:
                output = save_path + name + "_" + matrix + "_duration.svg"
                fig.savefig(output)
            if show:
                plt.show()


def plot_speed_up(duration, duration_seq, show=True, save=True, save_path=""):
    """Plots speed up charts for each algorithm for each matrix.

    Args:
        duration (dict): Data shaped like the output of
        `fill_parallel_duration`
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if show or save:
        for execution in duration:
            name, _, matrix = execution
            threads = duration[execution]['threads']
            time_transpose = duration[execution]['transpose']
            time_transpose_tr = duration[execution]['transpose_tr']
            time_ref = duration_seq[matrix]
            for time in time_transpose:
                time = time_ref / time
            for time in time_transpose_tr:
                time = time_ref / time
            fig, axes = plt.subplots(figsize=(8, 8))
            # , xlim=(threads[0], threads[-1])
            axes.set(xlabel='threads',
                     xticks=range(threads[0], threads[-1] + 1),
                     ylabel='duration (s)', title=name + " on " + matrix)
            axes.grid(True, linestyle='-.')
            line_transpose, = axes.plot(
                threads, time_transpose, linestyle=(0, (5, 7)), marker='o')
            line_transpose_tr, = axes.plot(threads, time_transpose_tr, 'ro--')
            line_transpose.set_label("transpose")
            line_transpose_tr.set_label("transpose_tr")
            axes.legend()
            if save:
                output = save_path + name + "_" + matrix + "_speed_up.svg"
                fig.savefig(output)
            if show:
                plt.show()


def minima(x):
    """Returns the minimum of `x` and the set of argmin.

    Args:
        x : 

    Returns:
        (, set): The first variable is the minimum of `x`. The second variable
        is the set of argmin.
    """
    minimum = x[0]
    indices = set()
    for i in range(len(x)):
        if x[i] < minimum:
            minimum = x(i)
            indices = set()
            indices.add(i+1)
        elif x[i] == minimum:
            indices.add(i+1)
    return (minimum, indices)


def find_minima(duration):
    """Finds minima (according to the number of threads) of each execution of
    algorithm in `duration` and returns a dictionary where the keys are the
    executions and values are the minima.

    Args:
        duration (dict): Data shaped like the output of
        `fill_parallel_duration`

    Returns:
        dict: Keys are execution of algorithm, values are the best duration and
        the according number of threads.
    """
    best = dict()
    for execution in duration:
        time_transpose = duration[execution]['transpose']
        time_transpose_tr = duration[execution]['transpose_tr']
        value = []
        time, threads = minima(time_transpose)
        value.append({'transpose': time, 'threads': threads})
        time, threads = minima(time_transpose_tr)
        value.append({'transpose_tr': time, 'threads': threads})
        best[execution] = value
    return best


def keep_faster_parallel(data_parallel, best):
    """Returns a copy of `data_parallel` where only the best execution
    (according to the number of threads) of each algorithm is kept.

    Args:
        data_parallel (list of dict): Data shaped like the output of `split`.
        best (dict): Data shaped like the output of `find_minima`.

    Returns:
        list of dict: A copy of `data_parallel` where only the best execution
    (according to the number of threads) of each algorithm is kept.
    """
    faster_parallel = list()
    for row in data_parallel:
        key = (row['name'], row['cflags'], row['matrix'])
        if row['threads'] in best[key][0]['threads']:
            faster_parallel.append(row)
    return faster_parallel


def plot_boxplot(data, available, show=True, save=True, save_path=""):
    """Plots boxplot for each algorithm for each matrix.

    Args:
        data (list of dict): Data shaped like the output of `csv_read`.
        available (set): Data shaped like the output of `csv_read`.
        show (bool, optional): Whether or not showing plots. Defaults to True.
        save (bool, optional): Whether or not saving plots. Defaults to True.
        save_path (str, optional): The path where charts are saved. Defaults to
        "".
    """
    if show or save:
        matrixs = set()
        for (_, _, matrix, _) in available:
            matrixs.add(matrix)
        for matrix in matrixs:
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.set(xlabel='name',
                   ylabel='duration (s)',
                   title=matrix)
            ax.grid(True, linestyle='-.')
            time_transpose = list()
            time_transpose_tr = list()
            labels = list()
            N_repeat = data[0]['total_step']
            j = 0
            for i in range(0, len(data), N_repeat):
                if data[i]['matrix'] == matrix:
                    time_transpose.append([0]*N_repeat)
                    time_transpose_tr.append([0]*N_repeat)
                    labels.append(data[i]['name'])
                    for k in range(N_repeat):
                        time_transpose[j][k] = data[i+k]['transpose']
                        time_transpose_tr[j][k] = data[i+k]['transpose_tr']
                    j = j + 1
            # print(time_transpose)
            bplot = ax.boxplot(time_transpose, sym="k+",  # patch_artist=True,
                               medianprops={'color': "blue"})
            bplot_tr = ax.boxplot(time_transpose_tr, sym="k+")
            bplot['medians'][0].set_label("transpose")
            bplot_tr['medians'][0].set_label("transpose_tr")
            plt.setp(ax, xticklabels=2*labels)
            ax.legend()
            # for patch in bplot_tr['boxes']:
            #     patch.set_facecolor("grey")
            if save:
                output = save_path + matrix + "_boxplot.svg"
                fig.savefig(output)
            if show:
                plt.show()


def create_ref_duration(data):
    """Returns a dictionary where keys are reference sequential algorithms in
    `data` and values mean durations.

    Args:
        data (list of dict): Data shaped like the output of `csv_read`

    Returns:
        dict: A dictionary where keys are reference sequential algorithms in
        `data` and values mean durations.
    """
    data_means = compute_means(data)  # list of dict
    duration = dict()
    for row in data_means:
        key = row["matrix"]  # (row['name'], row['cflags'], row["matrix"])
        value = {'compress': row['compress'], 'transpose': row['transpose'],
                 'compress_tr': row['compress_tr'],
                 'transpose_tr': row['transpose_tr']}
        duration[key] = value
    return duration


def main():
    filename = "benchmarks.csv"
    save_path = "charts/"
    argc = len(sys.argv)
    if argc == 2:
        filename = sys.argv[1]
        start_index = filename.find("/") + 1
        save_path = save_path + filename[start_index:-4] + "/"
        os.mkdir(save_path)
    elif argc > 2:
        print(
            "Usage: create_charts [input_filename].\nThe default filename is",
            filename)
        return
    data, available = csv_read(filename)
    data_sequential, data_parallel = split(data)
    # data_ref = csv_read("csv/classical_O2.csv")
    # duration_seq = create_ref_duration(data_ref)
    # data_parallel_means = compute_means(data_parallel)
    # duration = create_parallel_duration(available)
    # fill_parallel_duration(data_parallel_means, duration)
    # plot_duration(duration, save=False, save_path=save_path)
    # plot_speed_up(duration, duration_seq, save=False, save_path=save_path)
    # best = find_minima(duration)
    # faster_parallel = keep_faster_parallel(data_parallel, best)
    new_data = data  # data_sequential + faster_parallel
    plot_boxplot(new_data, available, show=False,
                 save=True, save_path=save_path)


if __name__ == "__main__":
    main()
