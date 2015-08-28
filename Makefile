default:
	gcc -c -Werror -fpic -Wall src/bcql.c -o bcql.o -lopenblas
	gcc -shared -o libbcql.so bcql.o -lopenblas

install:
	#TODO use command "install" in place of "mkdir", "cp", "chmod", ecc..
	mkdir -p /usr/include/bcql
	cp src/bcql.h /usr/include/bcql/
	cp libbcql.so /usr/lib64/
	chmod 755 /usr/include/bcql/
	chmod 644 /usr/include/bcql/bcql.h
	chmod 755 /usr/lib64/libbcql.so
clean:
	echo "TODO"


