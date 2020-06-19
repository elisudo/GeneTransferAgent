'''
There are 278 protein fastas, where no GTA homologs were found.
This code extracts these fastas along with 278 TRUE GTA fastas to obtain equally sized groups
and copies these fastas into a separate folder for subsequent analysis
'''


import os
import shutil
import random



order = 'Rhizobiales'
Busco_threshold = 80
write = True


list_true_and_false_GTA_accessions = []
list_GTA_found = []
list_No_GTA_found = []
list_all_accessions = []
list_nothing_found = []
list_unknown_GTA = []
BUSCO_dict = {}

os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\BUSCO_output')
with open('BUSCO_summary.txt','r') as f:
    BUSCO_lines = f.readlines()
    BUSCO_lines.pop(0)
    for BUSCO_line in BUSCO_lines:
        words = BUSCO_line.split('\t')
        BUSCO_dict[words[0]] = float(words[1])

for i in os.listdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\GTA_results'):
    if i.startswith('GCA'):
        accession = i[0:15]
        list_true_and_false_GTA_accessions.append(accession)
        
for i in os.listdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Prot_fastas'):
    if i.startswith('GCA'):
        accession = i[0:15]
        list_all_accessions.append(accession)

os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\GTA_results')
with open('GTA_clustering_results.csv', 'r') as f:
    GTA_output_lines = f.readlines()

for acc in list_true_and_false_GTA_accessions:
    for line in GTA_output_lines:
        words = line.split(',')
        if acc in words[0]:
            if 'True' in words[2]:
                list_GTA_found.append(acc)
            if 'False' in words[2]:
                list_No_GTA_found.append(acc)
            if 'unknown' in words[2]:
                list_unknown_GTA.append(acc)

for acc2 in list_all_accessions:
    if acc2 not in list_GTA_found and acc2 not in list_No_GTA_found and acc2 not in list_unknown_GTA:
        list_nothing_found.append(acc2)



genomes_without_GTA = []
for i in list_No_GTA_found:
    try:
        if BUSCO_dict[i] >= Busco_threshold:
            genomes_without_GTA.append(i)
    except KeyError:
        print('KeyError')
for i in list_nothing_found:
    try: 
        if BUSCO_dict[i] >= Busco_threshold:
            genomes_without_GTA.append(i)
    except KeyError:
        print('KeyError')

genomes_with_GTA = []
for i in list_GTA_found:
    try:
        if BUSCO_dict[i] >= Busco_threshold:
            genomes_with_GTA.append(i)
    except KeyError:
        print('KeyError')


if write == True:
    for i in genomes_without_GTA:
        original1 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Prot_fastas\\{i}.faa'
        target1 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Prot_fastas\\Subset_without_GTA\\{i}.faa'
        original2 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Nuc_fastas\\{i}.fna'
        target2 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Nuc_fastas\\Subset_without_GTA\\{i}.fna'
        shutil.copyfile(original1,target1)
        shutil.copyfile(original2,target2)
    
    GTA_containing_subset = random.sample(genomes_with_GTA,len(genomes_without_GTA))
    for i in GTA_containing_subset:
        original1 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Prot_fastas\\{i}.faa'
        target1 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Prot_fastas\\Subset_with_GTA\\{i}.faa'
        original2 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Nuc_fastas\\{i}.fna'
        target2 = f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Nuc_fastas\\Subset_with_GTA\\{i}.fna'
        shutil.copyfile(original1,target1)
        shutil.copyfile(original2,target2)