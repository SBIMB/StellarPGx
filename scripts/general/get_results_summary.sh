#!/bin/bash


while getopts s:o: option
do
case "${option}"
in

s) FILE1=${OPTARG};;
o) FILE2=${OPTARG};;

esac
done


for i in $(cat $FILE1); do 

    if [[ $(grep "Possible" ${i}*alleles) && $(grep "*" ${i}*alleles) ]]; then
    paste <(echo ${i}) <(grep "*" ${i}*alleles) <(grep "Possible" ${i}*alleles) | column -s $'\t' -t -o $'\t' >> $FILE2


    elif [[ $(grep "*" ${i}*alleles) ]]; then
    paste <(echo ${i}) <(grep "*" ${i}*alleles) | column -s $'\t' -t -o $'\t' >> $FILE2
    

    elif [[ $(grep "Possible" ${i}*alleles) ]]; then
    paste <(echo ${i}) <(grep "Possible" ${i}*alleles) | column -s $'\t' -t -o $'\t' >> $FILE2

    else
    paste <(echo ${i}) <(echo No_call) | column -s $'\t' -t -o $'\t' >> $FILE2
    fi

done
