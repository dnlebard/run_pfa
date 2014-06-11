#!/usr/bin/env python

import os
import subprocess
import fnmatch
import shutil
import sys

################# ARGUMENT PARSING ##############################
EXPECTED_ARGS = 12

num_args = len(sys.argv)
if num_args < EXPECTED_ARGS:
    print '\nUsage: %s pfa_exec pfa_inp_file pfa_output_file num_replicas num_proc num_frames data_dir output_dir pfa_inp_dir psf_file nfo_file' % sys.argv[0]
    print '\nEXAMPLE: %s /Users/dnlebard/code/bin/pfa.cuda.e pfa.inp pfa.out 20 4 200 /Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/data /Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/output/replicas /Users/dnlebard/cg-peptide-project/pep250/pfa-test pep250.psf pep250.nfo\n' % sys.argv[0]
    sys.exit()

pfa_exec   = sys.argv[1]
pfa_input  = sys.argv[2]
pfa_output = sys.argv[3]
num_replicas = int(sys.argv[4])
num_proc = int(sys.argv[5])
num_frames = sys.argv[6]
data_dir =  sys.argv[7] # OLD VERSION: '/Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/data'
output_base = sys.argv[8] # OLD VERSION:  '/Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/output/replicas'
pfa_inp_dir = sys.argv[9] # (contains psf file)... OLD VERSION: '/Users/dnlebard/cg-peptide-project/pep250/pfa-test'
psf_file_name = sys.argv[10]
nfo_file_name = sys.argv[11]
#################################################################


# Cluster analysis input (SHOULD THE CUTOFF BE TAKEN FROM THE COMMAND LINE?)
cluster_beads = 'GB GB1 GB2 GB3 AB\n'
cluster_res = 'AMY\n'
cluster_cmd = 'PSF CLUSTER GYR\n'
cluster_cutoff = '10.0\n'
cluster_file_out = 'cluster.dat'

for rep_idx in range(num_replicas):
    print '\n\n=========== REPLICA %d =====================' % rep_idx
    print '\nLOOPING OVER REPLICA #%d' % rep_idx

    # Make file selection based on the replica id
    file_by_rep_idx = '*rep' + str(rep_idx) + '*.dcd'

    # make a dcd file list
    dcd_file_list = []
    for root, dirs, files in os.walk(data_dir):
        for file_name in fnmatch.filter(files, file_by_rep_idx):
            dcd_file_list.append(os.path.join(root, file_name))

    # sort dcd files so they are in acending order of the ns in the path
    dcd_file_list.sort(key=lambda dcd_file_path: int(dcd_file_path.split('/')[-2].split('ns')[0]))

    # create the file name for the output dcd file
    last_dcd = dcd_file_list[-1]
    last_ns = dcd_file_list[-1].split('/')[-2]
    output_dcd = 'rep' + str(rep_idx) + '.' + last_ns + '.dcd'
    #print 'THE CONSTRUCTED OUTPUT dcd FILE: %s\n\n' % output_dcd


    str_dcd_file_list = '  '.join(dcd_file_list)
    #print '%s' % ' '.join(dcd_file_list)

    # go to output_dir
    output_dir = output_base + '/' + str(rep_idx)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)

    # Run the catdcd command
    catdcd_cmd = '~/bin/catdcd -o ' + output_dcd + ' ' + str_dcd_file_list
    print 'RUN CATDCD COMMAND .... '
    print 'catdcd command: %s ' % catdcd_cmd
    os.system(catdcd_cmd)
    print '....DONE!\n'

    # Create the nfo file
    nfo_file = open(nfo_file_name,"w")
    nfo_file.write('1\n')
    nfo_file.write('1\n')
    nfo_file.write('1\n')
    nfo_file.write('1\n')
    nfo_file.close()
    
    # Create psf.conf input file
    psf_conf_file = open('psf.conf',"w")
    psf_conf_file.write('0      0\n')
    psf_conf_file.write('0      0\n')
    psf_conf_file.write('0      0\n')
    psf_conf_file.write('0      0\n')
    psf_conf_file.write('1      5000\n')
    psf_conf_file.write('5000\n')
    psf_conf_file.close()

    # Create box.in input file
    box_inp_file = open('box.in',"w")
    box_inp_file.write('90.0000000\n')
    box_inp_file.write('90.0000000\n')
    box_inp_file.write('90.0000000\n')
    box_inp_file.write('209.254601\n')
    box_inp_file.write('209.254601\n')
    box_inp_file.write('209.254601\n')
    box_inp_file.close()

    # copy psf file
    shutil.copy(os.path.join(pfa_inp_dir, psf_file_name),output_dir)

    # Create PFA input file
    pfa_inp_file = open(pfa_input,"w")
    pfa_inp_file.write('box.in\n')  # Box file
    pfa_inp_file.write('ORTHO\n')   # ortho shape
    pfa_inp_file.write(cluster_cmd) # PFA COMMANDS
    pfa_inp_file.write((nfo_file_name + '\n')) # nfo file
    pfa_inp_file.write((psf_file_name  + '\n'))    # psf file name
    pfa_inp_file.write((output_dcd + '\n')) # dcd file name
    pfa_inp_file.write('DCD\n') # DCD file type
    pfa_inp_file.write('NO\n') # NO MPI IO
    pfa_inp_file.write(num_frames) # Number of frames in the file
    pfa_inp_file.write(cluster_beads) # CG BEAD NAMES
    pfa_inp_file.write(cluster_res) # RESIDUE NAME
    pfa_inp_file.write(cluster_cutoff)  # Cutoff length
    pfa_inp_file.write((cluster_file_out + '\n')) # Cluster output file name
    pfa_inp_file.write('CPU\n')  # CPU-only calculation
    pfa_inp_file.close()

    # Make and execute PFA command
    pfa_cmd = 'mpiexec -np %d %s < %s' % (num_proc,pfa_exec,pfa_input)
    print 'RUN PFA COMMAND .... '
    print 'pfa command: %s ' % pfa_cmd
    print '....DONE!'

    #subprocess.check_call('dos2unix mgo_input', stdout=open('/dev/null','w'), shell=True)
    #shutil.copy2(mgo_input, 'INPUT')
    #subprocess.check_call('mpirun -np 8 Pcrystal',

    # Remove dcd file
    os.remove(output_dcd)
    print '\n========================================\n'
