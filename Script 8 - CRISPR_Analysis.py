'''
1. go through all folders that start with GCA
2. in each folder, extract required info
3. make summary file
4. make overall summary file which includes all accession numbers and sum of all spacers in genomes with a complete system
'''

# importing the neccessary stuff
import os

## initializing some values and user specified arguments
order = 'Genomes_of_interest'
with_or_without = 'with and without' # Do you want to analyze the subset with or without the GTAs?

if order == 'MEGARHO':
    if with_or_without == 'with':
        folder = '\\Subset_with_GTA'
    elif with_or_without == 'without':
            folder = '\\Subset_without_GTA'
    else:
        print("please select either 'with' or 'without'")
else:
    folder = ''
## openening folder and initializing output file
os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\CRISPRCasFinder_results{folder}')
with open('TOTAL.CRISPR.summary.tsv', 'w') as f:
    f.write('Assembly.Accession\tCRISPR.Evidence_level\tCas.genes\tSystem\tNb.of.unique.spacers\tComplete.CRISPR\n')
f.close()

os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\CRISPRCasFinder_results{folder}')
with open('TOTAL.CAS.summary.csv', 'w') as f4:
    f4.write('Assembly.Accession,CAS-type\n')
f4.close()
    
# Going through all the files
total_nb_of_complete_systems = 0
total_nb_of_cas = 0
total_nb_of_files = 0
for i in os.listdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\CRISPRCasFinder_results{folder}'):
    
    TypeIV = False
    TypeIU = False
    TypeIID = False
    TypeIIIB = False
    TypeIIIA = False
    TypeIIC = False
    TypeIF = False
    TypeIE = False
    TypeIC = False
    TypeIA = False
                
    if i.startswith("GCA_"):
        accession = i[0:15]
        total_nb_of_files += 1
        os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\CRISPRCasFinder_results{folder}\\{i}\\TSV')
        try:
            with open('Crisprs_REPORT.xls','r') as f1:
                lines1 = f1.readlines()
        except FileNotFoundError:
            print(f'file not found for {i}, continuing to the next one')
            continue
        list_of_spacer_numbers = []
        list_of_duplicate_spacers = []
        total_spacer_count = 0
        total_duplicate_spacer_count = 0
        for line1 in lines1:
            if len(line1) <= 5:
                continue
            else:
                words1 = line1.split('\t')
                list_of_spacer_numbers.append(words1[14])
                list_of_duplicate_spacers.append(words1[3])
        if len(list_of_spacer_numbers) > 1:
            list_of_spacer_numbers2 = list_of_spacer_numbers[1:]
            for spacer_count in list_of_spacer_numbers2:
                total_spacer_count += int(spacer_count)
        if len(list_of_duplicate_spacers) > 1:
            list_of_duplicate_spacers2 = list_of_duplicate_spacers[1:]
            for duplicate_spacer_count in list_of_duplicate_spacers2:
                total_duplicate_spacer_count += int(duplicate_spacer_count)
       
        with open('CRISPR-Cas_summary.tsv','r') as f2:
            lines2 = f2.readlines()
        
        os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\CRISPRCasFinder_results{folder}')
        line_counter2 = 0
        with open ('TOTAL.CAS.summary.csv', 'a') as f4:
            Cas_list = []
            for line2 in lines2:
                words3 = line2.split('\t')
                try:
                    if 'CAS-Type' in line2 and line_counter2 != 0:
                        Cas_string = words3[4]
                        Cas_string =Cas_string.split(',')
                        Cas_string.pop()
                        for Cas in Cas_string: 
                            Cas2 = Cas.split('[')
                            Cas_list.append(Cas2[0])
                            f4.write(f'{i},{Cas2[0]}\n')
                            total_nb_of_cas += 1
#                            print('Cas detected')
                except IndexError:
                    continue
                except TypeError:
                    continue
                except ValueError:
                    continue
                finally:
                    line_counter2 += 1
        
        with open ('CRISPR.summary.csv', 'w') as f3:
           f3.write('Sequence_ID,Evidence_level,cas_number,Nb_of_spacers,cas_type\n')
        f3.close()
        
        Has_cas = False
        is_complete = False
        highest_level = 0
        line_counter = 0
        with open ('CRISPR.summary.csv', 'a') as f3:
            for line2 in lines2:
                if 'Type' in line2 and line_counter != 0:
                    Has_cas = True
                words2 = line2.split('\t')
                try:
                    if int(words2[2]) >= 1:
                        evidence_levels_and_text = words2[1]
                        evidence_levels1 = evidence_levels_and_text.split(',')
                        evidence_levels1.pop()
                        evidence_levels2 = []
                        for level in evidence_levels1:
                            evidence_levels2.append(int(level[-2]))
                        evidence_levels2.sort()
                        highest_level = evidence_levels2[-1]                      
                        f3.write(f'{words2[0]},{highest_level},{words2[5]},{total_spacer_count},{Cas_list}\n')
#                        print('write successful')
                except IndexError:
                    continue
                except TypeError:
                    continue
                except ValueError:
                    continue
                finally:
                    line_counter += 1
                f1.close()
                f2.close()
                f3.close()
            if highest_level == 4 and Has_cas == True:
                is_complete = True
                total_nb_of_complete_systems += 1
            
        
        if 'CAS-TypeIV' in Cas_list:
            TypeIV = True
        if 'CAS-TypeIU' in Cas_list:
            TypeIU = True
        if 'CAS-TypeIIID' in Cas_list:
            TypeIID = True
        if 'CAS-TypeIIIB' in Cas_list:
            TypeIIIB = True
        if 'CAS-TypeIIIA' in Cas_list:
            TypeIIIA = True
        if 'CAS-TypeIIC' in Cas_list:
            TypeIIC = True
        if 'CAS-TypeIF' in Cas_list:
            TypeIF = True
        if 'CAS-TypeIE' in Cas_list:
            TypeIE = True
        if 'CAS-TypeIC' in Cas_list:
            TypeIC = True
        if 'CAS-TypeIA' in Cas_list:
            TypeIA = True
            
        
        
            # write them to a new file which is the big ass summary
        os.chdir(f'C:\\Users\\elisc\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\eli\\BEP\\{order}\\CRISPRCasFinder_results{folder}')
        with open('TOTAL.CRISPR.summary.tsv', 'a') as f:
            f.write(f'{accession}\t{highest_level}\t{Has_cas}\t{Cas_list}\t{total_spacer_count - total_duplicate_spacer_count}\t{is_complete}\n')
            f.close()
print('\n\n\n---- REPORT ----')        
print(f'\n{total_nb_of_files} files were processed for {order} in the subset {with_or_without} GTAs')
print(f'{total_nb_of_complete_systems} complete CRISPR systems were detected ({total_nb_of_complete_systems / total_nb_of_files * 100:.2f}%)')
print(f'{total_nb_of_cas} Cas systems were detected')



