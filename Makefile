default:
	gcc -c -Werror -fpic -Wall src/libquat.c -o libquat.o -lopenblas
	gcc -shared -o libquat.so libquat.o

install:
	mkdir /usr/include/libquat
	cp src/libquat.h /usr/include/libquat/
	cp libquat.so /usr/lib64/

clean:
	echo "TODO"


