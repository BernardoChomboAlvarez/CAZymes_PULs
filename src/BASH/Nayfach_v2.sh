cd ~/THEEND/data

LIST=$(cat 1069MAGS_list.txt)

#Creation of each MAG directory and each MAG diamond analysis
for i in ${LIST[@]};do 
	mkdir $i
	mv fna/$i.fna.gz* $i/.

	#Execution of diamond.sh for each MAG
	gzip -d $i/*.gz | ln -s ~/scripts/diamond.sh $i/. | qsub -V $i/diamond.sh &

done

#Extraction of PUL ID and its discrmination based on obtained e-value
for i in ${LIST[@]};do
	#e-value <= 1e-4 and sort usign only PUL's core ID
	cut -f2,11 $i/pul.matches.tsv | tr '_' '\t' | cut -f1,3 | tr '\t' ';' | awk -F ";" '$2~/^[1-9]\.[0-9]+[e]/ {print}' | tr '-' '\t' | awk -F "\t" '{ if ($2 >= 4) print}' | tr '\t' '-' | tr ';' '\t' | sort | uniq -c | sort | uniq -c >> $i/$i.matches_evalue.out
	cut -f2,11 $i/pul.matches.tsv | tr '_' '\t' | cut -f1,3 | tr '\t' ';' | awk -F ";" '$2~/^[1-9]\.[0-9]+[e]/ {print}' | tr '-' '\t' | awk -F "\t" '{ if ($2 >= 4) print}' | tr '\t' '-' | tr ';' '\t' | cut -f1 | sort | uniq -c | awk '{print $2,$1}' >> $i/$i.matches_counts.out

done
