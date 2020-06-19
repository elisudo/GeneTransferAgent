cd /home/eli/GTA-Hunter-v1-master/; 
for genome in /home/eli/BEP/Rhodobacterales/Prot_fastas/*faa;
do python /home/eli/GTA-Hunter-v1-master/GTA_Hunter.py -g data/training/gta/3_gta.faa -v data/training/viral/3_viral.faa \
	-w data/training/gta/3_gta.dist data/training/viral/3_viral.dist \
	-q $genome \
	-o /home/eli/BEP/Rhodobacterales/Prot_fastas/GTA_results/ \
	--kmer --blast --min;
	mv /home/eli/BEP/Rhodobacterales/Prot_fastas/GTA_results/results.out $genome.out;
	mv $genome.out /home/eli/BEP/Rhodobacterales/Prot_fastas/GTA_results/;
	rm /home/eli/BEP/Rhodobacterales/Prot_fastas/GTA_results/blast.out;
	rm /home/eli/BEP/Rhodobacterales/Prot_fastas/GTA_results/*.faa;
	rm /home/eli/BEP/Rhodobacterales/Prot_fastas/GTA_results/no_homologs.txt;
done;