#!/bin/sh 

NAME=`/usr/bin/id -un`;
TMP="/tmp/cgiknot-$NAME";

/usr/bin/find $TMP/ -type f -a \! -mmin -60 -exec /bin/rm -f \{\} \;
