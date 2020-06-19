cd /data/eli/Subset_with_GTA/;
for accession in *.fna;
do cd /data/eli/CRISPRCasFinder/;
	perl CRISPRCasFinder.pl -in /data/eli/Genomes_of_interest/$accession \
	-out ./RESULTS4/$accession \
	-levelMin 0 \
	-cas ;
done;