#!/usr/bin/env python3

import csv
import itertools
import matplotlib.pyplot as plt
import numpy as np
import sys

plt.rcParams.update({'font.size': 12})

save_path = "charts/"
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
            converted_row = {'name': row['name'],
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
            data.append(converted_row)
            key = (converted_row['name'], converted_row['cflags'],
                   converted_row['matrix'], converted_row['threads'])
            available.add(key)
    return data, available


def compute_means(data):
    N_repeat = int(data[0]['total_step'])
    means = list()
    for i in range(0, len(data), N_repeat):
        duration = [0, 0, 0, 0]
        for j in range(0, N_repeat):
            duration[0] += data[i+j]['compress']
            duration[1] += data[i+j]['transpose']
            duration[2] += data[i+j]['compress_tr']
            duration[3] += data[i+j]['transpose_tr']
        for j in range(0, 4):
            duration[j] = duration[j] / N_repeat
        row = {'name': data[i]['name'], 'cflags': data[i]['cflags'],
               'cxxflags': data[i]['cxxflags'], 'threads': data[i]['threads'],
               'matrix': data[i]['matrix'], 'compress': duration[0],
               'transpose': duration[1], 'compress_tr': duration[2],
               'transpose_tr': duration[3]
               }
        means.append(row)
    return means


def is_sequential(name, cflags):
    is_gustavson = name == "Gustavson"
    is_std_sort = name == "std::sort"
    is_mkl_sequential = (name[:3] == "MKL") and (name[4:] == "sequential")
    is_sequential = is_gustavson or is_std_sort or is_mkl_sequential
    return is_sequential


def split(data):
    sequential = list()
    parallel = list()
    for row in data:
        if is_sequential(row['name'], row['cflags']):
            sequential.append(row)
        else:
            parallel.append(row)
    return sequential, parallel


def create_parallel_duration(available):
    duration = dict()
    # value = {'threads': list(), 'compress': list(), 'transpose': list(),
    #                  'compress_tr': list(), 'transpose_tr': list()}
    value = {'threads': list(), 'transpose': list(), 'transpose_tr': list()}
    for key in available:
        name, cflags, matrix, _ = key
        short_key = (name, cflags, matrix)
        if not is_sequential(name, cflags):
            duration[short_key] = value
    return duration


def fill_parallel_duration(data_parallel_means, duration):
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


def plot_speed_up(duration, show=True, save=True):
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
            output = save_path + name+"_"+matrix+"_speed_up.svg"
            if save:
                fig.savefig(output)
            if show:
                plt.show()


def minima(x):
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
    faster_parallel = list()
    for row in data_parallel:
        key = (row['name'], row['cflags'], row['matrix'])
        if row['threads'] in best[key][0]['threads']:
            faster_parallel.append(row)
    return faster_parallel


def plot_boxplot(data, available, show=True, save=True):
    if show or save:
        matrixs = set()
        for (name, cflags, matrix, threads) in available:
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
            bplot = ax.boxplot(time_transpose, sym="k+", #patch_artist=True,
                               medianprops={'color': "blue"})
            bplot_tr = ax.boxplot(time_transpose_tr, sym="k+")
            bplot['medians'][0].set_label("transpose")
            bplot_tr['medians'][0].set_label("transpose_tr")
            plt.setp(ax, xticklabels=2*labels)
            ax.legend()
            # for patch in bplot_tr['boxes']:
            #     patch.set_facecolor("grey")
            output = save_path + matrix + "_boxplot.svg"
            if save:
                fig.savefig(output)
            if show:
                plt.show()


def main():
    filename = "benchmarks.csv"
    argc = len(sys.argv)
    if argc == 2:
        filename = sys.argv[1]
    elif argc > 2:
        print("Usage: create_charts [input_filename].\nThe default filename is", filename)
        return
    data, new_available = csv_read(filename)
    # print(data)
    # print(available)
    # data_sequential, data_parallel, available, all_available = split(data)
    data_sequential, data_parallel = split(data)
    # for row in data_sequential:
    #     print("seq", row)
    # print("para", data_parallel)
    # # print('all',all_available)
    # data_parallel_means = compute_means(data_parallel)
    # for row in data_parallel_means:
    #     print(row)
    # print(data_parallel_means)
    # duration = create_parallel_duration(new_available)
    # fill_parallel_duration(data_parallel_means, duration)
    # plot_speed_up(duration, save=False)
    # best = find_minima(duration)
    # print(best)
    # faster_parallel = keep_faster_parallel(data_parallel, best)
    # for row in faster_parallel:
    #     print(row)
    # for row in itertools.chain(data_sequential, faster_parallel):
    #     print(row)
    new_data = data # data_sequential + faster_parallel
    plot_boxplot(new_data, new_available, show=False, save=True)


if __name__ == "__main__":
    main()
