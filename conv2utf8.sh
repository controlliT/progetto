#!/bin/bash

#Faccio un backup
cp progetto_3A_utf8.m progetto_3A_utf8.m.bak
cp progetto_3A.m progetto_3A.m.bak
cp progetto_3A_new_utf8.m progetto_3A_new_utf8.m.bak
cp progetto_3A_new.m progetto_3A_new.m.bak

echo "Creato backup"

#Conversione
iconv -f "ISO-8859-14" -t "UTF-8" progetto_3A_utf8.m -o progetto_3A.m
iconv -f "ISO-8859-14" -t "UTF-8" progetto_3A_new_utf8.m -o progetto_3A_new.m

echo "conversione completata"
