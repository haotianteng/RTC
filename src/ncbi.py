#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 23:15:44 2018

@author: haotianteng
"""

from Bio import Entrez
from ena import taxon
import thread


    
def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid.  we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace("_", "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList']

def Entrez_init(email = "havens.teng@gamil.com"):
    Entrez.email=email

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)

def iter_get_taxid(full_sp_n):
    full_name = full_sp_n
    full_sp_n = full_sp_n.strip()
    full_sp_n = full_sp_n.replace('_','+').replace(' ','+')
    tax_id = get_tax_id(full_sp_n)
    while (len(tax_id)==0) and len(full_sp_n.split('+'))>1:
        full_sp_n = '+'.join(full_sp_n.split('+')[:-1])
        tax_id = get_tax_id(full_sp_n)
    if len(tax_id) >3:
        print(full_name) 
    return tax_id

def get_domain(species_name):
    ###Iteratively search the species and return the suspect domain name
    tax_id= iter_get_taxid(species_name)
    tax_data=[get_tax_data(x) for x in tax_id ]
    hits = list()
    for idx, data in enumerate(tax_data):
        name = list()
        name.append(tax_data[idx][0]['ScientificName'])
        current_domain = data[0]['Lineage'].split(';')[1].strip()
        if ('OtherNames' in data[0].keys()) and ('Synonym' in data[0]['OtherNames']):
            name.append(data[0]['OtherNames']['Synonym'])
        hits.append({'taxid':tax_id[idx],'name':name,'domain':current_domain})
    return hits

    
def get_lineage(species_name):
    ###Iteratively search the species and return the suspect domain name
    tax_id= iter_get_taxid(species_name)
    tax_data=[get_tax_data(x) for x in tax_id ]
    hits = list()
    for idx, data in enumerate(tax_data):
        name = list()
        name.append(tax_data[idx][0]['ScientificName'])
        current_lineage = data[0]['Lineage'].strip().split(';')
        if ('OtherNames' in data[0].keys()) and ('Synonym' in data[0]['OtherNames']):
            name+=data[0]['OtherNames']['Synonym']
        hits.append({'taxid':tax_id[idx],'name':name,'lineage':current_lineage})
    return hits
