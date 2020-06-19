import os
import itertools
import operator
import statistics as st
import matplotlib.pyplot as plt

### THINGS YOU SHOULD CHANGE ###
order = 'Rhizobiales'
max_dif = 3000
min_GTA_genes = 10

### THINGS YOU SHOULDN'T CHANGE ###
# Initialixing some values and the main output file
GTA_genes = ['2','3','4','5','6','8','9','12','13','14','15'] # These are the key GTA genes that should be in a cluster
final_dict = {}
Sequence_ID_dict = {}
total_no_files = 0
total_clusters_found = 0

if order == 'Rhodobacterales':
    os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\pAGOS')
    with open('All defense systems.csv', 'r') as f2:
        lines2 = f2.readlines()
        f2.close()
    
os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\GTA_results')
with open('GTA_clustering_results.csv', 'w') as f:
    f.write('Assembly.Accession,no.of.homologs,GTA.Cluster,Sequence ID,Start_coord,Stop_coord,Length,Has_pAGOS,pAGOS_GTA_dist\n')
    f.close()

pAGOS_distances = []
CRISPR_distances = []
BREX_distances = []
Hachiman_distances = []
Thoeris_distances = []
ABI_distances = []
TA_distances = []
Shedu_distances = []
Septu_distances = []
Druantia_distances = []
RM_distances = []
Gabija_distances = []
Zorya_distances = []
Disarm_distances = []
Wadjet_distances = []
Solo_cas_4_distances = []
CBASS_distances = []

# Going through all the files
for i in os.listdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\GTA_results'):
    if i.endswith(".faa.out"):
        os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\GTA_results')
        total_no_files += 1
        accession = i[0:15] # Reading the accession from the title
        gene_numbers = []
        GTA_homologs = []
        Sequence_IDs = []
        gene_start = []
        gene_stop = []
        file = open(i,'r') # open file in read mode
        lines = file.readlines()
        file.close()
        line_index = 0
        while line_index < len(lines):
            if lines[line_index][0] == '>' and 'GTA' in lines[line_index]: # check whether a line contains a gene classified as GTA homolog
                # obtain sequence ID
                Sequence_IDs.append(lines[line_index][1:11])
                # obtain gene number
                gene_numbers.append(lines[line_index][12:16])
                # obtain GTA homolog
                if lines[line_index][-3] == ' ':
                    GTA_homologs.append(lines[line_index][-2:-1])
                else:
                    GTA_homologs.append(lines[line_index][-3:-1])
                # obtain gene_start and gene_stop
                line = lines[line_index]
                words = line.split()
                gene_start.append(int(words[2]))
                gene_stop.append(int(words[4]))
                # to avoid counting double GTA homologs, skip the next line in file
                line_index += 2
            else:
                line_index += 1
                
        # assuming that there is no overlap, the gene indices will be the same in both lists after sorting
        Sequence_IDs = [x for _,x in sorted(zip(gene_start,Sequence_IDs))]
        gene_start.sort()
        gene_stop.sort()
        distances = []
        for i in range(len(gene_start)-1):
            distance = abs(gene_start[i+1] - gene_stop[i]) # the distance between the end of the first gene and the start of the second gene
            distances.append(distance)
        
        cons = []
        for i in range(len(distances)):
            if distances[i] < max_dif:
                cons.append(i)
        
        # detecting clusters
        all_clusters = []
        for a, b in itertools.groupby(enumerate(cons), lambda i_x: i_x[0]-i_x[1]):
            all_clusters.append((list(map(operator.itemgetter(1), b))))
        
        # checking whether found clusters meet thresholds as set by user
        for i in range(len(all_clusters)):
            cluster_size = len(all_clusters[i]) + 1
            if cluster_size >= min_GTA_genes:
                judge2 = True
                total_clusters_found += 1
                start_coord = min(gene_start[all_clusters[i][0]], gene_start[all_clusters[i][-1] + 1])
                stop_coord = max(gene_stop[all_clusters[i][0]], gene_stop[all_clusters[i][-1] + 1])
                Sequence_ID = Sequence_IDs[all_clusters[i][0]]
                print(f'{accession} has RcGTA-like cluster; The cluster size is {cluster_size} genes, {abs(stop_coord-start_coord)}')
                break
            elif cluster_size >= 2 and cluster_size <=9:
                print (f'{accession} has {cluster_size} out of 11 genes in the cluster. Therefore, no call can be made')
                judge2 = 'unknown'
            else:
                print(f'{accession} has {cluster_size} out of 11 genes in the cluster. Therefore, there is no GTA cluster present')
                judge2 = False
            final_dict[accession] = judge2
        
        if order == 'Rhodobacterales':
            file = open(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\Nuc_fastas\\{accession}.fna')
            text = file.read()
            Sequences = text.split('>')
            for j in range(len(Sequences)):
                length = 0
                lines = Sequences[j].split('\n')
                ID = Sequences[j][0:10]
                for line in range(1,len(lines)):
                    length += len(lines[line])
                Sequence_ID_dict[ID] = length
            file.close()
            
        
            # Checking distance to pAGOS
            Has_pAGOS = False
            pAGOS_GTA_dist = 0
            Has_CRISPR = False
            CRISPR_GTA_dist = 0
            Has_BREX = False
            BREX_GTA_dist = 0
            Has_Thoeris = False
            Thoeris_GTA_dist = 0
            Has_Hachiman = False
            Hachiman_GTA_dist = 0
            Has_ABI = False
            ABI_GTA_dist = 0
            Has_TA = False
            TA_GTA_dist = 0
            Has_Septu = False
            Septu_GTA_dist = 0
            Has_Shedu = False
            Shedu_GTA_dist = 0
            Has_Druantia = False
            Druantia_GTA_dist = 0
            Has_RM = False
            RM_GTA_dist = 0
            Has_Gabija = False
            Gabija_GTA_dist = 0
            Has_BREX = False
            BREX_GTA_dist = 0
            Has_Zorya = False
            Zorya_GTA_dist = 0
            Has_Disarm = False
            Disarm_GTA_dist = 0
            Has_Wadjet = False
            Wadjet_GTA_dist = 0
            Has_Solo_cas_4 = False
            Solo_cas_4_GTA_dist = 0
            Has_CBASS = False
            CBASS_GTA_dist = 0
            
            for line in lines2:
                words = line.split(',')
                if accession in words[0] and words[3] == 'Pagos':
                    Has_pAGOS = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        pAGOS_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        pAGOS_distances.append(pAGOS_GTA_dist)
                elif accession in words[0] and 'CRISPR' in words[3]:
                    Has_CRISPR = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        CRISPR_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        CRISPR_distances.append(CRISPR_GTA_dist)
                elif accession in words[0] and 'BREX' in words[3]:
                    Has_BREX = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        BREX_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        BREX_distances.append(BREX_GTA_dist)
                elif accession in words[0] and 'Hachiman' in words[3]:
                    Has_Hachiman = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Hachiman_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Hachiman_distances.append(Hachiman_GTA_dist)
                elif accession in words[0] and 'Thoeris' in words[3]:
                    Has_Thoeris = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Thoeris_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Thoeris_distances.append(Thoeris_GTA_dist)
                elif accession in words[0] and 'ABI' in words[3]:
                    Has_ABI = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        ABI_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        ABI_distances.append(ABI_GTA_dist)
                elif accession in words[0] and 'Toxin-antitoxin' in words[3]:
                    Has_TA = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        TA_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        TA_distances.append(TA_GTA_dist)
                elif accession in words[0] and 'Shedu' in words[3]:
                    Has_Shedu = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Shedu_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Shedu_distances.append(Shedu_GTA_dist)
                elif accession in words[0] and 'Septu' in words[3]:
                    Has_Septu = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Septu_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Septu_distances.append(Septu_GTA_dist)
                elif accession in words[0] and 'Druantia' in words[3]:
                    Has_Druantia = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Druantia_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Druantia_distances.append(Druantia_GTA_dist)
                elif accession in words[0] and 'RM' in words[3]:
                    Has_RM = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        RM_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        RM_distances.append(RM_GTA_dist)
                elif accession in words[0] and 'Gabija' in words[3]:
                    Has_Gabija = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Gabija_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Gabija_distances.append(Gabija_GTA_dist)
                elif accession in words[0] and 'Zorya' in words[3]:
                    Has_Zorya = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Zorya_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Zorya_distances.append(Zorya_GTA_dist)
                elif accession in words[0] and 'DISARM' in words[3]:
                    Has_Disarm = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Disarm_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Disarm_distances.append(Disarm_GTA_dist)
                elif accession in words[0] and 'Wadjet' in words[3]:
                    Has_Wadjet = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Wadjet_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Wadjet_distances.append(Wadjet_GTA_dist)
                elif accession in words[0] and 'Solo cas4' in words[3]:
                    Has_Solo_cas_4 = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        Solo_cas_4_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        Solo_cas_4_distances.append(Solo_cas_4_GTA_dist)
                elif accession in words[0] and 'CBASS' in words[3]:
                    Has_CBASS = True
                    if Sequence_ID in words[1]:
                        dist1 = abs(int(words[5]) - start_coord)
                        dist2 = abs(int(words[6]) - start_coord)
                        dist3 = abs(int(words[5]) - stop_coord)
                        dist4 = abs(int(words[6]) - stop_coord)
                        CBASS_GTA_dist = min(dist1, dist2, dist3, dist4)/Sequence_ID_dict[Sequence_ID] * 100
                        CBASS_distances.append(CBASS_GTA_dist)
            
        
        
      
        
        # write output to file
        if order == 'Rhodobacterales':
            with open('GTA_clustering_results.csv', 'a') as f:
                f.write(f'{accession},{cluster_size},{judge2},{Sequence_ID},{start_coord},{stop_coord},{abs(stop_coord-start_coord)},{Has_pAGOS},{pAGOS_GTA_dist}\n')
                f.close()
        else:
            with open('GTA_clustering_results.csv', 'a') as f:
                f.write(f'{accession},{cluster_size},{judge2}\n')
                f.close()

print(f'\n\n---- REPORT ----\n{total_no_files} were analysed and {total_clusters_found} GTA clusters were discovered')
if order == 'Rhodobacterales':
    print(f'pAGOS:\n\tAverage distance:{sum(pAGOS_distances)/len(pAGOS_distances):.0f}\n\tSTDEV:{st.stdev(pAGOS_distances):.0f}')
    print(f'CRISPR:\n\tAverage distance:{sum(CRISPR_distances)/len(CRISPR_distances):.0f}\n\tSTDEV:{st.stdev(CRISPR_distances):.0f}')
    print(f'BREX:\n\tAverage distance:{sum(BREX_distances)/len(BREX_distances):.0f}\n\tSTDEV:{st.stdev(BREX_distances):.0f}')
    print(f'Hachiman:\n\tAverage distance:{sum(Hachiman_distances)/len(Hachiman_distances):.0f}\n\tSTDEV:{st.stdev(Hachiman_distances):.0f}')
    #print(f'Thoeris:\n\tAverage distance:{sum(Thoeris_distances)/len(Thoeris_distances):.0f}\n\tSTDEV:{st.stdev(Thoeris_distances):.0f}')
    print(f'ABI:\n\tAverage distance:{sum(ABI_distances)/len(ABI_distances):.0f}\n\tSTDEV:{st.stdev(ABI_distances):.0f}')
    print(f'TA:\n\tAverage distance:{sum(TA_distances)/len(TA_distances):.0f}\n\tSTDEV:{st.stdev(TA_distances):.0f}')
    print(f'Shedu:\n\tAverage distance:{sum(Shedu_distances)/len(Shedu_distances):.0f}\n\tSTDEV:{st.stdev(Shedu_distances):.0f}')
    print(f'Septu:\n\tAverage distance:{sum(Septu_distances)/len(Septu_distances):.0f}\n\tSTDEV:{st.stdev(Septu_distances):.0f}')
    print(f'Druantia:\n\tAverage distance:{sum(Druantia_distances)/len(Druantia_distances):.0f}\n\tSTDEV:{st.stdev(Druantia_distances):.0f}')
    print(f'RM:\n\tAverage distance:{sum(RM_distances)/len(RM_distances):.0f}\n\tSTDEV:{st.stdev(RM_distances):.0f}')
    print(f'Gabija:\n\tAverage distance:{sum(Gabija_distances)/len(Gabija_distances):.0f}\n\tSTDEV:{st.stdev(Gabija_distances):.0f}')
    print(f'Zorya:\n\tAverage distance:{sum(Zorya_distances)/len(Zorya_distances):.0f}\n\tSTDEV:{st.stdev(Zorya_distances):.0f}')
    #print(f'Disarm:\n\tAverage distance:{sum(Disarm_distances)/len(Disarm_distances):.0f}\n\tSTDEV:{st.stdev(Disarm_distances):.0f}')
    #print(f'Wadjet:\n\tAverage distance:{sum(Wadjet_distances)/len(Wadjet_distances):.0f}\n\tSTDEV:{st.stdev(Wadjet_distances):.0f}')
    print(f'Solo cas4:\n\tAverage distance:{sum(Solo_cas_4_distances)/len(Solo_cas_4_distances):.0f}\n\tSTDEV:{st.stdev(Solo_cas_4_distances):.0f}')
    print(f'CBASS:\n\tAverage distance:{sum(CBASS_distances)/len(CBASS_distances):.0f}\n\tSTDEV:{st.stdev(CBASS_distances):.0f}')
    
    data = [pAGOS_distances, CRISPR_distances, BREX_distances, Hachiman_distances, Thoeris_distances, ABI_distances, TA_distances, Shedu_distances, Septu_distances, Druantia_distances, RM_distances, Gabija_distances, Zorya_distances, Disarm_distances, Wadjet_distances, Solo_cas_4_distances, CBASS_distances]
    fig1, ax1 = plt.subplots()
    ax1.set_title('Prokaryotic defense systems distance to GTA cluster')
    ax1.boxplot(data)
    plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], ['pAGOS', 'CRISPR', 'BREX','Hachiman','Thoeris','ABI','TA','Shedu','Septu','Druantia','RM','Gabija','Zorya','DISARM','Wadjet','Solo cas4','CBASS'],rotation=90)
    plt.ylim((0,100))
    plt.ylabel('Distance expressed as % of chromosome length')

