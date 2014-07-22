#!/usr/bin/python
import csv, os, glob
import pickle
from BeautifulSoup import BeautifulSoup
 
 
"""
# Shell script for fetching data from CAZy database 
#!/bin/bash
 
for i in $(echo {A..W})
do
wget -r -l 1 --random-wait -A b[0-9]*.html http://www.cazy.org/b$i.html
done

# Adopted and modified based on works pulled from https://gist.github.com/jstvz/1386797
"""

cazy_dict = {}

# Set working directory
WorkingDir = ""
RawHTMLDir = "/raw_html"

os.chdir(WorkingDir + RawHTMLDir)
for filename in glob.iglob('*.html'):
    print filename
    fam_list = []
    fh = open(filename, 'rb')
    fh_out_str = WorkingDir + filename.strip('.html') + '.xls'
    fh_out = open(fh_out_str, 'wb')
    soup = BeautifulSoup(fh)
    
    # Get the organism name, taxid, and lineage. Initialize the org container
    org_info = soup.find('div', {'class' : 'org'})
    
    try:
        org_name = org_info.find('font', {'class' : 'titre_cazome'}).renderContents()
        taxid = org_info.find('a', {'target' : 'ncbitaxid'}).renderContents()
        lineage_string = str(org_info).split(":")[-1].split('<')[0].strip()
        fam_list = [('Organism', org_name), ("NCBI_TaxID", taxid), ('Org_Lin',
            lineage_string)]
    except AttributeError:
        continue
    
    # Parse the tables, skipping the final protein list
    t = soup.findAll('table')
    for table in t[0:-1]:
        head = table.findAll('th')
        family_name = head[0].findAll('div',{'class':"titre_famille"})[0].renderContents()
        rows = table.findAll('tr')
        for tr in rows:
            cols = tr.findAll('td')
            for col in cols:
                # """
                # descend into each table column and try to pick out the data by
                # CSS class.
                # E.g: '<td class="classe" id="navigationtab2"><a
                # href="http://www.cazy.org/CE1.html" class="nav2"><div
                # class="famille">1</div></a><div class="nombre_de">2</div></td>'
                # """
                family = col.find("div", {"class" : "famille"}).renderContents()
                count = col.find("div", {"class" : "nombre_de"}).renderContents()
                link_out = col.find("a",{"class" : "nav2"}, href=True)['href']
                family_abbrv = link_out.split('/')[-1].strip('.html')
                fam_list.append((family_abbrv , count))
    # Write the dict to a new csv filename
    cazy_dict[org_name] = dict(fam_list)
    first_keys = ['Organism', 'NCBI_TaxID', 'Org_Lin'] # cazy_dict[org_name].keys()
    last_keys = [k for k in cazy_dict[org_name].keys() if k not in first_keys]
    last_keys.sort()
    [first_keys.append(k) for k in last_keys]
    dict_writer = csv.DictWriter(fh_out, first_keys, dialect = 'excel-tab',
            extrasaction = 'ignore')
    dict_writer.writer.writerow(first_keys)
    dict_writer.writerow(cazy_dict[org_name])
    fh_out.close()
 
 
 
with open('/media/haldane/trab/data/cazy/data/output_cazy.pkl', 'wb') as cazy_dict_out:
    pickle.dump(cazy_dict, cazy_dict_out)
    
    


