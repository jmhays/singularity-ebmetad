CPPFLAGS=-I./ 
SRC=*.cpp
OUT=libplumed.so

default: $(OUT)

links:
	./do_links.sh

$(OUT): links
	g++  $(CPPFLAGS) -O3 -DPLUMED_GENERIC_CXX -fPIC --shared -o $(OUT) *.cpp

clean:
	rm -f $(OUT) 
	find -type l -delete
	find -name \*~ -delete

veryclean: clean
	rm *~ TAGS
