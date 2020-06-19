
import os
from Bio import Entrez
import time
from termcolor import colored

order = 'Genomes_of_interest'
Entrez.email = 'eli.schalkers@gmail.com'

os.chdir('C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\MEGARHO\\Anvio')
# file can be found online at https://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt
with open('assembly_summary_genbank.txt', 'r', encoding="utf8") as f:
    TAXID_lines = f.readlines()

os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\MEGARHO\\GTA_results')
with open('GTA_clustering_results.csv', 'r') as f:
    GTA_lines_temp = f.readlines()
    GTA_lines1 = []
    for line in GTA_lines_temp:
        GTA_lines1.append(line[0:-1])
        
os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\Rhizobiales\\GTA_results')
with open('GTA_clustering_results.csv', 'r') as f:
    GTA_lines_temp = f.readlines()
    GTA_lines2 = []
    for line in GTA_lines_temp:
        GTA_lines2.append(line[0:-1])        

GTA_lines = GTA_lines1 + GTA_lines2
    
with open(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Anvio\\namedict.txt', 'w') as f:
    f.write('name\tcontigs_db_path\n')
    for fasta in os.listdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Nuc_fastas'):
        if fasta.startswith('GCA'):
            f.write(f'{fasta[0:15]}\t{fasta[0:15]}.db\n')
#    for fasta in os.listdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Prot_fastas\\Subset_with_GTA'):
#        if fasta.startswith('GCA'):
#            f.write(f'{fasta[0:15]}\t{fasta[0:15]}.db\n')
            
os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Anvio')
with open('namedict.txt', 'r') as f:
    anvio_lines = f.readlines()
    anvio_lines.pop(0)
    
with open('Additional_info.txt','w') as f:
    f.write('genome_id\tOrder\tFamily\tGenus\tSpecies\tGTA\tNumber_of_ORFs\n')
    
accessions = []
GTA_logical = False
for anvio_line in anvio_lines:
    segments = anvio_line.split()
    accessions.append(segments[0])

taxid_dict = {}
order_dict = {}
family_dict = {}
genus_dict = {}
species_dict = {}
ORF_dict = {}

Index_error_counter1 = 0
Index_error_counter2 = 0
Index_error_counter3 = 0
Index_error_counter4 = 0

failed_accessions = []

for accession in accessions:
    with open(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Prot_fastas\\{accession}.faa') as f:
        data = f.read()
        ORF_dict[accession] = len(data.split('>'))
    for TAXID_line in TAXID_lines:
        if accession in TAXID_line:
            words = TAXID_line.split('\t')
            taxid_dict[accession] = words[6]
    try:
        handle = Entrez.efetch('taxonomy', id=taxid_dict[accession], rettype='xml')
        response = Entrez.read(handle)
    except KeyError:
        break
    for entry in response:
#        sci_name = entry.get('ScientificName')
        lineage_taxa = entry.get('Lineage').split(';')
        try:
            order_dict[accession] = lineage_taxa[4]
        except IndexError:    
            print(colored("No lineage data found",'red'))
            Index_error_counter1 += 1
            order_dict[accession] = None
            family_dict[accession] = None
            genus_dict[accession] = None
            species_dict[accession] = None
            continue
        try:
            family_dict[accession] = lineage_taxa[5]
        except IndexError:    
            print(colored("Family level not found, but we've got the rest so don't worry!",'red'))
            Index_error_counter2 += 1
            family_dict[accession] = None
            genus_dict[accession] = None
            species_dict[accession] = None
            continue
        try:    
            genus_dict[accession] = lineage_taxa[6]
        except IndexError:    
            print(colored("Genus level not found, but we've got the rest so don't worry!",'red'))
            Index_error_counter3 += 1
            genus_dict[accession] = None
            species_dict[accession] = None
            continue
        try:
            species_dict[accession] = lineage_taxa[7]
        except IndexError:    
            print(colored("Species level not found, but the rest we've got so don't worry!",'red'))
            Index_error_counter4 += 1
            species_dict[accession] = None
            continue
    for GTA_line in GTA_lines:
        words = GTA_line.split(',')
        if words[0] == accession:
            GTA_logical = words[2]
            break
        else:
            GTA_logical = False
    with open('Additional_info.txt','a') as f:
        try:
            f.write(f'{accession}\t{order_dict[accession]}\t{family_dict[accession]}\t{genus_dict[accession]}\t{species_dict[accession]}\t{GTA_logical}\t{ORF_dict[accession]}\n')
            print(f'accession: {accession}\torder: {order_dict[accession]}\tfamily: {family_dict[accession]}\tgenus: {genus_dict[accession]}\tspecies: {species_dict[accession]}\tGTA: {GTA_logical}\tNumber of ORFs: {ORF_dict[accession]}')
        except KeyError:    
            print('Key Error encountered')
            failed_accessions.append(accession)
    time.sleep(0.34)
print(f'Number of index_error1: {Index_error_counter1}')
print(f'Number of index_error2: {Index_error_counter2}')
print(f'Number of index_error3: {Index_error_counter3}')
print(f'Number of index_error4: {Index_error_counter4}')
