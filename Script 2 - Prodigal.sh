cd /home/eli/BEP/Rhodobacterales/Nuc_fastas/; 
for accession in *.fna; 
do prodigal -i /home/eli/BEP/Rhodobacterales/Nuc_fastas/$accession \
	-c -a /home/eli/BEP/Rhodobacterales/Prot_fastas/"$accession".faa \
	-o /home/eli/BEP/Rhodobacterales/Prodigal_annotations/$accession.gff \
	-f gff; 
done;