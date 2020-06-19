cd /home/eli/BEP/Rhodobacterales/;
for accession in $(cat /home/eli/BEP/Rhodobacterales/Assembly_Accessions.csv);
do esearch -db assembly -query $accession | elink -target nuccore -name \
        assembly_nuccore_insdc | elink -target protein | efetch -format fasta \
        > $accession.faa;
done;