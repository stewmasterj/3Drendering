#!/bin/bash
if [ -z $2 ]; then
 echo "usage:"
 echo "cropNcnvt.sh FILE.ppm SCR.rc"
 echo "  FILE.ppm     is dumped from View to convert and crop to PNG"
 echo "  SCR.rc       is the screen.rc file used in the rendering"
 exit
else
 dims=`grep subscreen $2 | awk '{print $5-$3+1"x"$4-$2+1"+"$3-1"+"$2-1}'`
 echo "using subscreen and crop string: $dims"
 convert -crop ${dims} $1 ${1%\.*}.png
fi
