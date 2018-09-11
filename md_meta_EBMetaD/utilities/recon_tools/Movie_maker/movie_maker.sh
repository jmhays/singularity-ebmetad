#!/bin/bash

# Information used to make movies
nrepeat=10              # Number of times to repeat clustering images
startframe=2476          # First frame in movie
endframe=2625           # Last frame in movie

# Information on trajectory
cluster_stride=25       # Number of frames between cluster analysis steps (this will depend on how often you choose to output frames from vmd)
tstride=1000            # Time between cluster analysis steps (in ps)
tstep=`echo $tstride \/ \( 1000*$cluster_stride \) | bc -l`  # Time between frames (in ns)

# Information on locations of the various images from which we will compose the 
# movie
BASDIR=/global/gtribello/Isoleucine/0.1kT/Movie/Basins         # directory in which to find basins
basinbase=basin                                                # Base for names of basin images 
TRAJDIR=/global/gtribello/Isoleucine/0.1kT/Movie/Traj_frames   # directory in which to find frames from trajectory
trajbase=isoleu                                                # Base for names of frame images

# Create the clustering label
convert -background white -fill black -font Helvetica-Bold -pointsize 48 -size 300x60 -gravity South label:CLUSTERING clustering.gif

nbasins=`wc -l BASINS | awk '{print $1}'`
nbasins=$(($nbasins-1))

# Look for frames directory and establish 
# how many frames are in the movie already
if [ -d Frames ] ; then
   frameno=`ls Frames/frame* | wc -l | awk '{print $1}'`
   franeno=$(($frameno-1))
   if [ -e Frames/frame.$frameno.jpg ] ; then
      echo Something is wrong file already exists
      exit
   fi
   checkf=$(($frameno-1))
   if [ ! -e Frames/frame.$checkf.jpg ] ; then
      echo Something is wrong - there are not enough files $frameno $checkf
      exit
   fi
else
   mkdir Frames
   frameno=0
fi

for i in `seq $startframe $endframe` ; do

  # Establish whether to trajectory analysis
  modcheck=`expr $i % $cluster_stride`
  if [ $i -eq 0 ] ; then
     modcheck=1
  fi
  echo Frame $i $modcheck 

  # Create the number of the file with the zeros in
  if [ $i -lt 10 ] ; then
    ii=000$i
  elif [ $i -lt 100 ] ; then
    ii=00$i
  elif [ $i -lt 1000 ] ; then
    ii=0$i
  else
    ii=$i
  fi
  
  # This deals with frames for which there is a trajectory analysis
  if [ $modcheck -eq 0 ] ; then
  
    # Establish when clusters analysis was done
    ncluster=`expr $i / $cluster_stride`
    ctime_u=$(($ncluster*$tstride))
    ctime_l=$(($ctime_u-$tstride))
  
    echo All basins generated between $ctime_l and $ctime_u
    echo -ne Basin numbers " "
    # Find basins in BASINS
    nbas=0
    for j in `seq 0 $nbasins` ; do
        # Create the number of the file with the zeros in
        if [ $j -lt 10 ] ; then
          jj=000$j
        elif [ $j -lt 100 ] ; then
          jj=00$j
        elif [ $j -lt 1000 ] ; then
          jj=0$j
        else 
          jj=$j
        fi
        ctime=`head -n $(($j+1)) BASINS | tail -1 | awk '{print $2}'`
        logcheck=`echo $ctime \<\= $ctime_u \&\& $ctime \> $ctime_l | bc -l`
        if [ $logcheck -eq 1 ] ; then
          # Use image magic to convert the basins to a jpg and resize
          convert $BASDIR/$basinbase.$jj.tga -resize '30%' -transparent "#FFFFFF" bas.$nbas.png 
          echo -ne $j " "
          nbas=$(($nbas+1))
        fi
        # This is not strictly necessary - just breaks out of the loop if we are 
        # past the final clustering time - makes this faster
        logcheck=`echo $ctime \> $ctime_u | bc -l`
        if [ $logcheck -eq 1 ] ; then
           break
        fi 
    done
    echo TOTAL $nbas

    # Create the label
    stime=`echo $i \* $tstep | bc -l`
    printf "time = %6.2f ns" $stime > label.dat
    # Create a picture of the timestep label
    convert -background white -fill black -font Helvetica-Bold -pointsize 40 label:@./label.dat label.gif

    # Create the labeled image
    composite -gravity South label.gif $TRAJDIR/$trajbase.$ii.tga label_frame.jpg 

    # This will assemble the image for the clustering step
    if [ $nbas -eq 0 ] ; then
       composite -gravity North clustering.gif label_frame.jpg clusterpic.jpg
    elif [ $nbas -eq 1 ] ; then
       composite -gravity North clustering.gif label_frame.jpg stage1.jpg
       composite -gravity NorthWest bas.0.png stage1.jpg clusterpic.jpg
    elif [ $nbas -eq 2 ] ; then
       composite -gravity North clustering.gif label_frame.jpg stage1.jpg
       composite -gravity NorthWest bas.0.png stage1.jpg stage2.jpg
       composite -gravity NorthEast bas.1.png stage2.jpg clusterpic.jpg
    elif [ $nbas -eq 3 ] ; then
       composite -gravity North clustering.gif label_frame.jpg stage1.jpg
       composite -gravity NorthWest bas.0.png stage1.jpg stage2.jpg
       composite -gravity NorthEast bas.1.png stage2.jpg stage3.jpg
       composite -gravity SouthWest bas.2.png stage3.jpg clusterpic.jpg   
    elif [ $nbas -eq 4 ] ; then 
       composite -gravity North clustering.gif label_frame.jpg stage1.jpg
       composite -gravity NorthWest bas.0.png stage1.jpg stage2.jpg
       composite -gravity NorthEast bas.1.png stage2.jpg stage3.jpg    
       composite -gravity SouthWest bas.2.png stage3.jpg stage4.jpg
       composite -gravity SouthEast bas.3.png stage4.jpg clusterpic.jpg
    else
      echo ERROR - too many basins found to put in a single frame - rewrite script
      exit 
    fi

    # Create the clustering frames for use in the movie
    for k in `seq 1 $nrepeat` ; do
        cp clusterpic.jpg Frames/frame.$frameno.jpg
        frameno=$(($frameno+1))
    done

    # Get rid of the crap ready for next creation
    rm -f clusterpic.jpg bas* stage* label_frame.jpg label.dat label.gif
  
  # This deals with other frames
  else
    # calculate the current time in nanoseconds
    stime=`echo $i \* $tstep | bc -l`
    printf "time = %6.2f ns" $stime > label.dat
    # Create a picture of the timestep label
    convert -background white -fill black -font Helvetica-Bold -pointsize 40 label:@./label.dat label.gif
    # Use image to add the label to the image
    composite -gravity South label.gif $TRAJDIR/$trajbase.$ii.tga Frames/frame.$frameno.jpg
    # convert $TRAJDIR/$trajbase.$ii.tga Frames/frame.$frameno.jpg 
    frameno=$(($frameno+1))
    # Remove stuff 
    rm -f label.dat label.gif
  fi

done

rm -f clustering.gif
