cd /home/eli/BEP/Rhodobacterales/Prot_fastas/
for genome in *.faa;
do cd /home/eli/BEP/Rhodobacterales/BUSCO_output;
	busco -i /home/eli/BEP/Rhodobacterales/Prot_fastas/$genome \
	-o $genome.folder \
	-m prot \
	-l rhodobacterales_odb10 \
	-c 10;
mv /home/eli/BEP/Rhodobacterales/BUSCO_output/${genome}.folder/run_rhodobacterales_odb10/short_summary.txt \
	/home/eli/BEP/Rhodobacterales/BUSCO_output/ ;
mv /home/eli/BEP/Rhodobacterales/BUSCO_output/short_summary.txt \
	/home/eli/BEP/Rhodobacterales/BUSCO_output/${genome}.txt;
rm -r /home/eli/BEP/Rhodobacterales/BUSCO_output/${genome}.folder;
done;


cd /home/eli/BEP/Rhodobacterales/BUSCO_output
for file in *.faa.txt;
do score=$(head -8 $file | tail -1 | tail -c+4 | head -c4);
	accession=$(echo $file | head -c15);
	printf "%s\t" $accession >> BUSCO_summary.txt;
	printf "%s\n" $score >> BUSCO_summary.txt
done;