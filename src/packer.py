#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 21:23:18 2018
This module will pack the fast5 reads into data batch according to the readmap.
@author: haotianteng
"""
import h5py,os,sys

def getID(fast5_f):
    with h5py.File(fast5_f,'r') as handle:
        protocal_id = handle['UniqueGlobalKey/tracking_id'].attrs['protocol_run_id']
    return protocal_id
    
def read_map(read_sp_map):
    """Read read->species dict from file. """
    read_sp = {}
    with open(read_sp_map ,'r') as h:
        for line in h:
            if line.startswith(">"):
                crt_sp = line.strip()[1:].strip(',')
            else:
                if crt_sp=="GRCh38":
                    continue
                read_sp[line.strip()] = crt_sp
    return read_sp

def read_kingdom(sp_kingdom_table):
    """Read species->kindom dict from file."""
    sp_kd = {}
    with open(sp_kingdom_table,'r') as h:
        for line in h:
            sp_line = line.strip().split(',')
            sp_kd[sp_line[0]] = sp_line[1]
    return sp_kd
    
def check_match(read_sp,sp_kd):
    """Check if read->species map accrospond with species->kingdom map"""
    miss_read = list()
    miss_sp = list()
    for read in read_sp.keys():
        if read_sp[read] not in sp_kd.keys():
            miss_read.append(read)
            if read_sp[read] not in miss_sp:
                miss_sp.append(read_sp[read])
    return miss_read,miss_sp

def fn_kingdom(fast5_folder,read_sp,sp_kingdom):
    """Build filename-> kingdom look-up table for all the fast5 files in the
    given directory and all child directory"""
    f_kd = dict()
    kd_f = dict()
    for cr_dir,_,file_list in os.walk(fast5_folder):
        for f in file_list:
            if f.endswith('fast5') or f.endswith('f5'):
                root = h5py.File(os.path.join(cr_dir,f),'r')
                read_id = list(root['/Raw/Reads'].values())[0].attrs['read_id']
                if read_id in read_sp.keys():
                    kd = sp_kingdom[read_sp[read_id]]
                f_kd[f] = kd
                if kd in kd_f.keys():
                    kd_f[kd].append(f)
                else:
                    kd_f[kd] = [f]
    return f_kd,kd_f
def run(args):
    if len(args)<5:
        print("Need at least 4 parameters: %s [read_species_map] [species_kingdom_f] [fast5_directory] [output_dir]")
    read_sp_f = args[1]
    sp_kingdom_f = args[2]
    f5_dir = args[3]
    o_dir = args[4]
    print("Read -> Species mapping.")
    read_sp = read_map(read_sp_f)
    print("Species -> kingdom mapping.")
    sp_kd = read_kingdom(sp_kingdom_f)
    print("Checking missing read and missing species in dict.")
    miss_read,miss_sp = check_match(read_sp,sp_kd)
    print("Following reads didn't found:")
    print(miss_read)
    print("Following species didn't found:")
    print(miss_sp)
    print("Filename -> kingdom.")
    f_kd,kd_f = fn_kingdom(f5_dir,read_sp,sp_kd)
    output_path =os.path.join(o_dir,'filename_kingdom.csv')
    print("Write output to %s"%(output_path))
    if not os.path.isdir(o_dir):
        os.mkdir(o_dir)
    with open(output_path,'w+') as o_h:
        for kd in kd_f.keys():
            o_h.write(">%s\n"%(kd))
            for fn in kd_f[kd]:
                o_h.write("%s\n"%(fn))
if __name__=="__main__":
#==============================================================================
#Test code
#==============================================================================
#    read_sp_f = '/media/Linux_ex/MetaGenomics_Data/species2reads.map'
#    sp_kingdom_f =  '/media/Linux_ex/MetaGenomics_Data/species_kingdom.csv'
#    print("Reading read-species mapping.")
#    read_sp = read_map(read_sp_f)
#    print("Reading species-kingdom mapping.")
#    sp_kd = read_kingdom(sp_kingdom_f)
#    print("Checking missing read and missing species in dict.")
#    miss_read,miss_sp = check_match(read_sp,sp_kd)
#    print("Following reads didn't found:")
#    print(miss_read)
#    print("Following species didn't found:")
#    print(miss_sp)
#    TEST_DIR = '/media/Linux_ex/MetaGenomics_Data/raw_data'
#    f_kd,kd_f = fn_kingdom(TEST_DIR,read_sp,sp_kd)

    run(sys.argv[:])
    
    