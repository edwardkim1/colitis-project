#!/bin/bash

currLoc="$PWD"
path="${currLoc}"
new_dir_name=$1

cd ${currLoc}

echo -n "extracting tar into ${new_dir_name}..."

#EXTRACT TAR IN DESIRED DIRECTORY AND GUNZIP EVERYTHING
mkdir $new_dir_name
tar -xf GSE134809_RAW.tar -C $path/$new_dir_name
cd $path/$new_dir_name
echo "done"
echo -n "gunzip..."
gunzip *
echo "done"
echo -n "organizing files..."
#MOVE FILES INTO 
for filename in *;
do
    extension="${filename##*.}"
    name=$(basename $filename)
    dir_name=${name%_*}
    if [[ ! -e $dir_name ]]; then
    	mkdir $dir_name
	mv ${dir_name}_genes.tsv $dir_name/genes.tsv
    	mv ${dir_name}_barcodes.tsv $dir_name/barcodes.tsv
    	mv ${dir_name}_matrix.mtx $dir_name/matrix.mtx
    fi
done
echo "done"