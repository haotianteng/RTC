#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 17:12:40 2018

@author: haotianteng
"""
from ena import taxon
from ncbi import get_lineage,Entrez_init
from threading import Thread
from Queue import Queue
from difflib import SequenceMatcher as SM
import time

Meta_path = "../test"
species_list = Meta_path + "/species_list.dat"
THREAD_NUM =6

def read_species_list(list_file):
    species_name = list()
    with open(list_file) as f:
        for line in f:
            species_name.append(line.strip()[1:].strip(','))
    return species_name

def search(i,q,l,f_l):
    """Search the strain name on the database and establish a species-domain
    look-up dictionary, first search on ENA database, if no exact match is 
    found, search on NCBI database, if still no exact match, the species and
    the search result with most similar name will be added into fail_list.
    """
    while q.qsize()>0:
        ex_match = False
        sp = q.get()
        print("[Thread:%d Queue:%d]Searching taxon information for %s on ENA database"%(i,q.qsize(),sp))
        ena_hits = taxon(sp)
        sm_score = list()
        #Sleep for 0.5 second to prevent the database from rejecting access.
        for idx,hit in enumerate(ena_hits):
            sp_norm = sp.replace('+',' ').replace('_',' ').replace('sp. ','')
            hit_norm = [name.replace('+',' ').replace('_',' ').replace('sp. ','') for name in hit['name']]
            idty_chck = any([sp_norm == name for name in hit_norm])
            if idty_chck:
                l[sp] = ena_hits[idx]
                ex_match = True
                q.task_done()
                print("[Thread:%d Queue:%d]Taxon information hasbee found for %s on ENA database"%(i,q.qsize(),sp))
                break
            else:
                sm_score.append(max([SM(None,sp_norm,name).ratio() for name in hit_norm]))
            
        if not ex_match:
            print("[Thread:%d Queue:%d]Exact match is not found in ENA database, searching %s on NCBI database"%(i,q.qsize(),sp))
            ncbi_hits = get_lineage(sp)
            for idx,hit in enumerate(ncbi_hits):
                hit_norm = [name.replace('+',' ').replace('_',' ').replace('sp. ','') for name in hit['name']]
                idty_chck = any([sp_norm == name for name in hit_norm])
                if idty_chck:
                    l[sp]=ncbi_hits[idx]
                    ex_match = True
                    q.task_done()
                    print("[Thread:%d Queue:%d]Taxon information hasbee found for %s on NCBI database"%(i,q.qsize(),sp))
                    break
                else:
                    sm_score.append(max([SM(None,sp_norm,name).ratio() for name in hit_norm]))
        if not ex_match:
            print("[Thread:%d Queue:%d]%s can't be found in both database, added into fail list."%(i,q.qsize(),sp))
            hits = ena_hits + ncbi_hits
            if len(hits)==0:
                f_l[sp]=[]
            else:
                print(sm_score)
                print(len(hits))
                max_idx = sm_score.index(max(sm_score))
                f_l[sp] = hits[max_idx]
            q.task_done()
        
if __name__=="__main__":
#==============================================================================
# Test code
#==============================================================================
#    test_species = "Burkholderia 383"
#    ncbi_hits = get_lineage(test_species)
#    ena_hits = taxon(test_species)
    
    
#==============================================================================
#Metagenomics species typing
#==============================================================================
    Entrez_init()
    hits_dict = dict()
    fail_dict = dict()
    sp_list = read_species_list(species_list)
    q = Queue()
    for sp in sp_list:
        q.put(sp)
    for i in range(THREAD_NUM):
        worker = Thread(target = search,args = (i,q,hits_dict,fail_dict))
        worker.setDaemon(True)
        worker.start()

    q.join()
#        