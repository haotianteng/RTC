#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 20:42:23 2018

@author: haotianteng
This module is build to help search on European Nucleotide Archive (ENA) database.
"""
import urllib2
import xml.etree.ElementTree as ET
def search(qstring,entrez,mode):
    response = urllib2.urlopen('http://www.ebi.ac.uk/ena/data/search?query=%s&result=%s&display=%s'%(qstring,entrez,mode))
    return response.read()

def taxon(strain_name):
    strain_name = strain_name.replace('_','+').replace(' ','+')
    xml_data = search(strain_name,'taxon','xml')
    root = ET.fromstring(xml_data)
    hits = list()
    for hit in root:
        name = list()
        name.append(hit.attrib['scientificName'])
        for child in hit:
            if child.tag=='synonym':
                name.append(child.attrib['name'])
        lineage = parse_lineage(hit[0])
        taxid = hit.attrib['taxId']
        hits.append({'name':name,'lineage':lineage,'taxid':taxid})
    return hits
def parse_lineage(lineage_element):
    """
    Input args:
        lineage_element: a xml Element by ET module.
    """
    lineage = dict()
    for taxon in lineage_element:
        if taxon.attrib['hidden'] =="true":
            continue
        try:
            lineage[taxon.attrib['rank']]=[taxon.attrib['scientificName'],taxon.attrib['taxId']]
        except KeyError:
#            print("Taxon entry %s does not has rank information, skipped."%(taxon.attrib['scientificName']))
            continue
    return lineage
    

if __name__=="__main__":    
    hits = taxon('Salinarchaeum_laminariae')
    