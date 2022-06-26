# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:35:00 2022

@author: Lucija
"""
import numpy as np

def get_three_bins(filename, bin_size_equal = False, limit1 = None, limit2 = None):
    try:    
        file = open("../sequences/seq_prediction/" + filename, 'r')
    except:
        return
    print("getting bins")
    print(filename)
    lines = file.readlines()
    line_segments = [line.replace(' ', '').replace('\n', '').split(',') for line in lines]
    sol_dict = [] 
    scaled_sol = [] 
    for line in line_segments:
        if line[0] == 'SEQUENCEPREDICTIONS':
            scaled_sol.append(float(line[3])) 
            sol_dict.append([
                line[1],
                float(line[2]),
                float(line[3]),
                float(line[4]),
                float(line[5])
            ]) 
    bin2_start = 1 / 3
    bin3_start = 2 / 3
    if bin_size_equal:
        bin2_start = np.quantile(scaled_sol, 1/3) 
        bin3_start = np.quantile(scaled_sol, 2/3)
    if limit1 != None:
        bin2_start = limit1
    if limit2 != None:
        bin3_start = limit2
    bins = [[],[],[]]
    for entry in sol_dict:
        if entry[2] < bin2_start:
            bins[0].append(entry[0])
            continue
        if entry[2] < bin3_start:
            bins[1].append(entry[0])
            continue
        bins[2].append(entry[0])
    print(bin2_start, bin3_start)
    print(len(bins[0]), len(bins[1]), len(bins[2]))
    file.close()
    return bins 

def get_three_bins_multiple_files(filenames,  bin_size_equal = False, limit1 = None, limit2 = None):
    lines = []
    print("getting multiple file bins")
    for filename in filenames:
        try:    
            file = open("../sequences/seq_prediction/" + filename, 'r')
        except:
            return
        print(filename)
        lines += file.readlines()
        file.close()
    line_segments = [line.replace(' ', '').replace('\n', '').split(',') for line in lines]
    sol_dict = [] 
    scaled_sol = []
    for line in line_segments:
        if line[0] == 'SEQUENCEPREDICTIONS':
            scaled_sol.append(float(line[3])) 
            sol_dict.append([
                line[1],
                float(line[2]),
                float(line[3]),
                float(line[4]),
                float(line[5])
            ])
    bin2_start = 1 / 3
    bin3_start = 2 / 3
    if bin_size_equal:
        bin2_start = np.quantile(scaled_sol, 1/3) 
        bin3_start = np.quantile(scaled_sol, 2/3)
    if limit1 != None:
        bin2_start = limit1
    if limit2 != None:
        bin3_start = limit2
    bins = [[],[],[]]
    for entry in sol_dict:
        if entry[2] < bin2_start:
            bins[0].append(entry[0])
            continue
        if entry[2] < bin3_start:
            bins[1].append(entry[0])
            continue
        bins[2].append(entry[0])
    print(bin2_start, bin3_start)
    print(len(bins[0]), len(bins[1]), len(bins[2]))
    return bins 
 
def write_bins(bins, filename_aligned, filename_unaligned):
    print("writing bins")
    print(filename_aligned, filename_unaligned)
    try:
        all_sequences_aligned = open("../sequences/training_validation/" + filename_aligned + ".fa", 'r')
    except:
        return
    lines_aligned = all_sequences_aligned.readlines()
    try:
        all_sequences_unaligned = open("../sequences/training_validation/" + filename_unaligned + ".fa", 'r')
    except:
        return
    lines_unaligned = all_sequences_unaligned.readlines()
    for bin_num in range(len(bins)):
        bin_write = bins[bin_num]
        file_aligned = open("../sequences/training_validation/" + filename_aligned + "_sol" + str(bin_num) + ".fa", 'w')
        file_unaligned = open("../sequences/training_validation/" + filename_unaligned + "_sol" + str(bin_num) + ".fa", 'w')
        data_write_aligned = ""
        data_write_unaligned = ""
        print(len(bin_write))
        for seq_id in bin_write:
            for line_num in range(0, len(lines_aligned), 2): 
                if lines_aligned[line_num].strip() == seq_id: 
                    data_write_aligned += lines_aligned[line_num].strip() + "\n" + lines_aligned[line_num + 1].strip() + "\n"
                    break
            for line_num in range(0, len(lines_unaligned), 2):
                if lines_unaligned[line_num].strip() == seq_id : 
                    data_write_unaligned += lines_unaligned[line_num].strip() + "\n" + lines_unaligned[line_num + 1].strip() + "\n"
                    break
        file_aligned.write(data_write_aligned)
        file_aligned.close()
        file_unaligned.write(data_write_unaligned)
        file_unaligned.close()
    all_sequences_aligned.close()
    all_sequences_unaligned.close()
    

def get_number_bins(filename, bin_size_equal = False, limit1 = None, limit2 = None):
    try:    
        file = open("../sequences/seq_prediction/" + filename, 'r')
    except:
        return
    print("getting and writing conditions")
    print(filename)
    lines = file.readlines()
    line_segments = [line.replace(' ', '').replace('\n', '').split(',') for line in lines]
    sol_dict = [] 
    scaled_sol = [] 
    for line in line_segments:
        if line[0] == 'SEQUENCEPREDICTIONS':
            scaled_sol.append(float(line[3])) 
            sol_dict.append([
                line[1],
                float(line[2]),
                float(line[3]),
                float(line[4]),
                float(line[5])
            ]) 
    bin2_start = 1 / 3
    bin3_start = 2 / 3
    if bin_size_equal:
        bin2_start = np.quantile(scaled_sol, 1/3) 
        bin3_start = np.quantile(scaled_sol, 2/3)
    if limit1 != None:
        bin2_start = limit1
    if limit2 != None:
        bin3_start = limit2
    string_bin = "" 
    for entry in sol_dict:
        if entry[2] < bin2_start: 
            string_bin += "0\n"
            continue
        if entry[2] < bin3_start:
            string_bin += "1\n" 
            continue
        string_bin += "2\n" 
    print(bin2_start, bin3_start) 
    print(string_bin.count("0"), string_bin.count("1"), string_bin.count("2"))
    file.close()
    file_conditions = open("../sequences/training_validation/" + filename.replace('seq_prediction', 'conditions'), 'w')
    file_conditions.write(string_bin)
    file_conditions.close() 
    
bins_val = get_three_bins('seq_prediction_ll_val.txt', limit1 = 0.34, limit2 = 0.447)
get_number_bins('seq_prediction_ll_val.txt', limit1 = 0.34, limit2 = 0.447)
write_bins(bins_val, 'luxafilt_llmsa_val', 'll_val')

bins_train = get_three_bins('seq_prediction_ll_train.txt', limit1 = 0.34, limit2 = 0.447)
get_number_bins('seq_prediction_ll_train.txt', limit1 = 0.34, limit2 = 0.447)
write_bins(bins_train, 'luxafilt_llmsa_train', 'll_train')

get_three_bins('seq_prediction_ll_val.txt')
get_three_bins('seq_prediction_ll_train.txt')
get_three_bins_multiple_files(['seq_prediction_ll_val.txt', 'seq_prediction_ll_train.txt'])

get_three_bins('seq_prediction_ll_val.txt', limit1 = 0.34, limit2 = 0.447)
get_three_bins('seq_prediction_ll_train.txt', limit1 = 0.34, limit2 = 0.447)
get_three_bins_multiple_files(['seq_prediction_ll_val.txt', 'seq_prediction_ll_train.txt'], limit1 = 0.34, limit2 = 0.447)
         
get_three_bins('seq_prediction_ll_val.txt', True)
get_three_bins('seq_prediction_ll_train.txt', True)
get_three_bins_multiple_files(['seq_prediction_ll_val.txt', 'seq_prediction_ll_train.txt'], True)