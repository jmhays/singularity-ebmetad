#!/usr/bin/awk -f 
#
BEGIN{tot_con=1}{ printf "%7s %4d %4d %4d %5.1f %2d %2d %5.1f %8.5f\n",$1,tot_con,$3,$4,$5,$6,$7,$8,$9; tot_con=tot_con+1}
