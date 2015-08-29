default:
	gcc -c -Werror -fpic -Wall src/bcsk.c -o bcsk.o -lopenblas
	gcc -shared -o libbcsk.so bcsk.o -lopenblas

install:
	#TODO use command "install" in place of "mkdir", "cp", "chmod", ecc..
	mkdir -p /usr/include/bcsk
	cp src/bcsk.h /usr/include/bcsk/
	cp libbcsk.so /usr/lib64/
	chmod 755 /usr/include/bcsk/
	chmod 644 /usr/include/bcsk/bcsk.h
	chmod 755 /usr/lib64/libbcsk.so
clean:
	echo "TODO"


