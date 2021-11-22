#!/bin/bash

chrlist="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

cd data/cadd
echo "Splitting CADD..."
for i in $chrlist;  
do   
awk '{
	if($1== "'"$i"'") 
		print $0; 
	else if($1> "'"$i"'") 
	exit
	}' hg19_cadd13.txt > cadd_$i.txt
done
python cadd_chunk.py -r $1
echo "CADD file split is complete."



cd ../dbnsfp
echo "Splitting dbnsfp..."

for i in $chrlist;  
do   
awk '{
	if($1== "'"$i"'") 
		print $0; 
	else if($1> "'"$i"'") 
	exit
	}' dbnsfp41a.txt > dbnsfp_$i.txt
done
python dbnsfp_chunk.py -r $1
echo "dbnsfp file split is complete."


cd ../revel
echo "Splitting revel..."
for i in $chrlist;  
do   
awk '{
	if($1== "'"$i"'") 
		print $0; 
	else if($1> "'"$i"'") 
	exit
	}' hg19_revel.txt > revel_$i.txt
done

python revel_chunk.py -r $1
echo "revel file split is complete."



cd ../mcap
echo "Splitting mcap..."
for i in $chrlist;  
do   
awk '{
	if($1== "'"$i"'") 
		print $0; 
	else if($1> "'"$i"'") 
	exit
	}' hg19_mcap13.txt > mcap_$i.txt
done
python mcap_chunk.py -r $1
echo "mcap file split is complete."
