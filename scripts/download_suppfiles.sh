#!/bin/bash

EMAIL=$2

while read -r suppfilename || [ -n "$suppfilename" ]
do
ID=$(echo ${suppfilename} | sed 's/^.*\(GSE[0-9]*\).*$/\1/g')
if [[ ! -f "suppl/${suppfilename}" ]]
then
curl -sS -H 'Expect:' \
-o "suppl/${suppfilename}" \
--user anonymous:${EMAIL} \
ftp://ftp.ncbi.nlm.nih.gov/geo/series/${ID:0:-3}nnn/${ID}/suppl/${suppfilename}
else 
echo "${suppfilename} is present"
fi
done < $1
