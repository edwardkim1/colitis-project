#!/bin/bash

currLoc="$PWD"
path="${currLoc}"
new_dir_name=$1

cd ${currLoc}

echo -n "extracting tar into ${new_dir_name}..."

#EXTRACT TAR IN DESIRED DIRECTORY AND GUNZIP EVERYTHING
mkdir $new_dir_name
tar -xf GSE125527_RAW.tar -C $path/$new_dir_name
cd $path/$new_dir_name
echo "done"
echo -n "gunzip..."
gunzip *
echo "done"
echo "extraction completed"
