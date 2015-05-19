#!/bin/bash

if [ -e $1 ]; then
	(head -n 1 $1 && tail -n +2 $1 | sort -rs -k3,3n -k2,2d -S 20%) > $1.tmp
	mv $1.tmp $1
fi
