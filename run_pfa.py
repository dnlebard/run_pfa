#!/usr/bin/env python

import os
import subprocess
import fnmatch
import shutil
import sys

################# ARGUMENT PARSING ##############################
EXPECTED_ARGS = 12

num_args = len(sys.argv)
if num_args != EXPECTED_ARGS:
    print '\nUsage: %s pfa_exec pfa_inp_file pfa_output_file num_replicas num_proc data_dir output_dir pfa_inp_dir psf_file nfo_file' % sys.argv[0]
    print '\nEXAMPLE: %s /Users/dnlebard/code/bin/pfa.cuda.e pfa.inp pfa.out 20 4 /Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/data /Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/output/replicas /Users/dnlebard/cg-peptide-project/pep250/pfa-test pep250.psf pep250.nfo num_beads\n' % sys.argv[0]
    sys.exit()

pfa_exec   = sys.argv[1]
pfa_input  = sys.argv[2]
pfa_output = sys.argv[3]
num_replicas = int(sys.argv[4])
num_proc = int(sys.argv[5])
data_dir =  sys.argv[6] # OLD VERSION: '/Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/data'
output_base = sys.argv[7] # OLD VERSION:  '/Users/dnlebard/cg-peptide-project/pep250/python_run_pfa/output/replicas'
pfa_inp_dir = sys.argv[8] # (contains psf file)... OLD VERSION: '/Users/dnlebard/cg-peptide-project/pep250/pfa-test'
psf_file_name = sys.argv[9]
nfo_file_name = sys.argv[10]
num_beads = sys.argv[11]
#################################################################


# Cluster analysis input (SHOULD THE CUTOFF BE TAKEN FROM THE COMMAND LINE?)
cluster_beads = 'GB GB1 GB2 GB3 AB\n'
cluster_res = 'AMY\n'
cluster_cmd = 'PSF CLUSTER GYR\n'
cluster_cutoff = '10.0\n'
cluster_file_out = 'cluster.dat'

#num_beads = 20000

for rep_idx in range(num_replicas):
    print '\n\n=========== REPLICA %d =====================' % rep_idx
    print '\nLOOPING OVER REPLICA #%d' % rep_idx

    # Make file selection based on the replica id
    file_by_rep_idx = '*rep' + str(rep_idx) + '.*.dcd'
    print 'THE VALUE OF file_by_rep_idx is "%s"' % file_by_rep_idx

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


    str_dcd_file_list = ' '.join(dcd_file_list)
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
    get_frames_cmd = '~/bin/catdcd ' + output_dcd + " |grep Total |awk '{print $3}'"
    #print 'GET FRAMES catdcd command: %s ' % get_frames_cmd
    num_frames = subprocess.check_output(get_frames_cmd, shell=True)
    #print 'NEW NUMBER OF FRAMES: %s' % num_frames

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
    psf_conf_file.write(('1      ' + num_beads + '\n'))
    psf_conf_file.write((num_beads + '\n'))
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
    #pfa_cmd = 'mpiexec -np %d %s < %s > %s' % (num_proc,pfa_exec,pfa_input,pfa_output)
    pfa_cmd = 'mpiexec -np %d %s' % (num_proc,pfa_exec)
    print 'RUNNING PFA COMMAND: |%s < %s >  %s| for replica %d of %d .... ' % (pfa_cmd,pfa_input,pfa_output,rep_idx,num_replicas)
    subprocess.call(pfa_cmd, stdout=open(pfa_output, 'w'), stdin=open(pfa_input), stderr=subprocess.STDOUT, shell=True)
    
    print '....DONE!'

    # Remove dcd file
    os.remove(output_dcd)
    print '\n========================================\n'
