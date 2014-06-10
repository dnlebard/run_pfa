#!/usr/bin/env python

import os
import fnmatch
data_dir = '/Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/data'
output_base = '/Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/output/replicas'

# TODO:
#       (1) Take the data dir from the command line args
#       (2) Take the output dir from the command line args


n_replicas = 20

for rep_idx in range(n_replicas):
    print '\nLOOPING OVER REPLICA #%d' % rep_idx

    # Make file selection based on the replica id
    file_by_rep_idx = '*rep' + str(rep_idx) + '*.dcd'

    #print 'file by idx : %s' % file_by_rep_idx

    # make a dcd file list
    dcd_file_list = []
    for root, dirs, files in os.walk(data_dir):
        for file_name in fnmatch.filter(files, file_by_rep_idx):
            dcd_file_list.append(os.path.join(root, file_name))


    #print 'THE LAST DCD FILE: %s\n\n' % dcd_file

    # sort dcd files so they are in acending order of the ns in the path
    dcd_file_list.sort(key=lambda dcd_file_path: int(dcd_file_path.split('/')[-2].split('ns')[0]))


    # create the file name for the output dcd file
    last_dcd = dcd_file_list[-1]
    last_ns = dcd_file_list[-1].split('/')[-2]
    output_dcd = 'rep' + str(rep_idx) + '.' + last_ns + '.dcd'
    print 'THE CONSTRUCTED OUTPUT dcd FILE: %s\n\n' % output_dcd


    str_dcd_file_list = ''.join(dcd_file_list)
    #print '%s' % ' '.join(dcd_file_list)

    catdcd_cmd = '~/bin/catdcd -o ' + output_dcd + ' ' + str_dcd_file_list

    #print 'RUN CATDCD COMMAND .... DONE!'
    print 'catdcd command: %s ' % catdcd_cmd

    output_dir = output_base + '/rep' + str(rep_idx)

    # Make output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
