#!/bin/bash
declare -A a
declare -A b
n=0
[[ -z $1 ]] && exit 1
[[ -z $2 ]] && exit 1
for t in $1$2.log
do 
	a["${n}"]="${t}"
	idx=$(echo "${t}"| sed -e 's@E\(.*\)_tol.*@\1@g')
	[[ $2 = "*" ]] && idx=$(echo "${t}"| sed -e 's@E.*_tol\(.*\)\.log@\1@g')
	b["${idx}"]=${n}
	((n++))
done

for t in $(echo "${!b[@]}" | sed -e "s@ @\n@g" | sort -g)
do
	n=${b["${t}"]}
	cat ${a["${n}"]}	
done
