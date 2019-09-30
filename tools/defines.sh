#!/bin/bash

grep -r "#if" src/* | cut -d " " -f 2- | cut -f1 -d"/" | sed -e 's/\ifdef//g' -e 's/\if//g' -e 's/\defined//g' -e 's/[()!#,|&\//<>=+*]//g' | awk -v RS=" " 'length()==0{next}{print $0}' | sed '/^[0-9]/ d' | tr '\r' '\n' | tr '\n\n' '\n' | sort -u > tmp_defines_used.txt
sed -e 's/\#//g' Template-Config.sh | grep -v "\-\-" | awk '{if ($1) print $1}' | awk -F"\=" '{print $1}' | sort > tmp_defines_exist.txt

comm -13 tmp_defines_exist.txt tmp_defines_used.txt

rm tmp_defines_used.txt
rm tmp_defines_exist.txt
