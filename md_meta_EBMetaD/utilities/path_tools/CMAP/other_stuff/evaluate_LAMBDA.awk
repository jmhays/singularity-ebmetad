#!/usr/bin/awk -f 
#
# Calculate LAMBDA 
#
{
if($NF!=0){ value[$1,$2]=$3 ;if($1>max){max=$1}}
}
END{ 
#
# calculate average distance i, i+1
#
 n=0
 mean=0.0
 mean2=0.0
 for(i=1;i<max;i++){
    for(j=i+1;j<=max;j++){
        if(j-1==i ){
          n++
          mean=mean + value[i,j]
          mean2=mean2 + value[i,j]**2
                   }
    }

 }
 mean=mean/n
 mean2=mean2/n
 stddev = ((mean2-mean**2)/n)**0.5
 printf(" \n")
 printf(" Average Distance i, i+1  ::  %f \n",mean)
 printf(" Standard Deviation       ::  %f \n",stddev)
 printf(" Percentage               ::  %f \n",(stddev/mean)*100)
 printf("** SUGGESTED LAMBDA       ::  %f \n",2.3/mean)
 printf(" \n ")
}

