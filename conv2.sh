#!/bin/bash

#Scrippettino per convertire i file in codifiche differenti. 

#Controllo argomenti
if [[ $# -ne 3 ]]; then
    echo "Parametri in ingresso non validi"
    echo "Usage: ./conv2ansi.sh <encoding from> <encoding to>"
    exit 1
fi

folder=$1
fromEncoding=$2
toEncoding=$3

if [[ ! -d "$folder" ]]; then
    echo "Il primo argomento deve essere una cartella"
    exit 1
fi

for $subFolder in "$folder"; do
    if [[ -d "$subFolder" ]]; then

        if [[ ! -d "$subFolder/backup" ]]; then
            mkdir "$subFolder/backup"
            echo "creata cartella backup in $subFolder"
        fi


        for $file in "$subFolder/*.m"; do
            #cp "$file" "$subFolder/backup/$file.bak"
            #iconv -f $fromEncoding -t $toEncoding "$file" -o "$file.$toEncoding.m"
        done
    fi
done