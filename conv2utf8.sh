#!/bin/bash

#Faccio un backup
cp progetto_3A_utf8.m progetto_3A_utf8.m.bak
cp progetto_3A.m progetto_3A.m.bak

echo "Creato backup"

#Conversione
iconv -f "windows-1252" -t "UTF-8" progetto_3A.m -o progetto_3A_utf8.m

echo "conversione completata"
