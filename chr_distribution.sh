#!usr/bin/env bash
#Author: Deepali L. Kundnani
#First positional argument: chr list
#Second positional argument: File location in String format
#

if [ "$1" = "" ]; then
    echo "Please enter the file/files with chrlist/chrsizes or sample bed file works"
    echo "Usage: bash chr_distribution.sh chr.sizes filelist out"
    exit 1
fi

if [ "$2" = "" ]; then
    echo "file with a list of files with their respective locations. Recommended to run it where the files ar epresent and give output in another folder if required"
    echo "Usage: bash chr_distribution.sh chr.sizes filelist out"
    exit 1
fi

if [ "$3" = "" ]; then
    echo "Default output file name: output"
    out='output'
    echo "Usage: bash chr_distribution.sh chr.sizes filelist out"
else
    out=$3
fi

if [ "$1" = "-h" ]; then
    echo "Usage: bash chr_distribution.sh chr.sizes filelist' out"
    exit 1
fi


#Count hte ribonucleotides in each chromosome
chrs=$(cut -f1 $1)
files=$(cat $2)
echo -e 'chr\tsize\t'$files >> $out
for chr in $chrs
do 
size=$(grep -e ${chr}'\s' $1)
arr=($size)
for file in $files
do
count=$(grep -e ${chr}'\s' *$file* | wc -l) 
arr+=($count)
done
echo ${arr[*]} >> $out
done

