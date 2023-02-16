#!/bin/bash
# Jeff Eilbott, 2017, jeilbott@surveybott.com

# inputs
FILE="${1}*"
FILE=$(echo $FILE | awk '{print $1}')
if [  -e "$FILE" ]; then
	ARGS=
	if [ -f "$FILE" ]; then
		ARGS="$(cat $FILE)"
	fi
	BASE=$(dirname $0)
	$BASE/ABA_bott.sh $ARGS
	rm $FILE
fi
