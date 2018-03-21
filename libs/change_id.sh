#1/bin/bash

cat $2 | \
awk -v OUT=$1 \
'BEGIN{
	id = 0
}
{
	if(NR%4==1){
		print "@"id > OUT
		id+=1
	}
	else{
		print $0 > OUT
	}
}
'