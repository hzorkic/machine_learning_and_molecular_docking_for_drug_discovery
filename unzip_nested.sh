# Unzips all files (within a folder) in a directory to specified location

for dir in */; do
    cd $dir
    for file in *gz; do
    	filename="$(cut -d'.' -f1 <<<$file)"
    	gunzip -c "${file}" > /stor/home/js88749/TCGA/RNAseq/"${filename}"
	done
    cd ..
done