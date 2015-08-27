default:
	gcc -c -Werror -fpic -Wall src/libquat.c -o libquat.o -lopenblas
	gcc -shared -o libquat.so libquat.o

install:
	mkdir -p /usr/include/libquat
	cp src/libquat.h /usr/include/libquat/
	cp libquat.so /usr/lib64/
	chmod 755 /usr/include/libquat/
	chmod 644 /usr/include/libquat/libquat.h
	chmod 755 /usr/lib64/libquat.so
clean:
	echo "TODO"


