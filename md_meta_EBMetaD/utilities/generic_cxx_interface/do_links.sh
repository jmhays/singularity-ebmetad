rm -f *.cpp *.h
for file in ../../common_files/*.c
do
  ln -s $file $(basename $file .c).cpp
done

ln -s ../../common_files/metadyn.h metadyn.h

