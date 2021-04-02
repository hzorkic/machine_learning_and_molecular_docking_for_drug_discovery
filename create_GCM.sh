# Extracts gene counts from all RNAseq files in a directory
i=1;
awk '{ print $1 }' $file > "temp$i">temp0
for file in *; do
	# skip over non rnaseq files 
	if [ "$file" == "temp" ] || [ "$file" = "COAD_GCM" ] ; then
    	continue;
	fi
	awk '{ print $2 }' $file > "temp$i"
    i=$((i+1))
done

# Combines all temp files to matrix
paste temp[0-55] | column -s $'\t' -t >COAD_GCM.csv