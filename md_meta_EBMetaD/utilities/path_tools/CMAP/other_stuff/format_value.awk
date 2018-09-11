#!/usr/bin/awk -f 
#
BEGIN{tot_con=1}{ if($1!="END") {printf "%4d %4d %4d    %12.8f\n",tot_con,$2,$3,$4; tot_con=tot_con+1} else {printf "%3s\n",$1; tot_con=1}}
