openfhe:
	mkdir -p /usr/local/include/openfhe/build && cd /usr/local/include/openfhe/build && cmake -DBUILD_UNITTESTS=OFF -DBUILD_EXAMPLES=OFF -DBUILD_BENCHMARKS=OFF .. && make -j && make install && ldconfig


# make LWE-KEM -B
LWE-KEM:
	cd /app/target ; rm LWE-KEM ; g++ -o LWE-KEM LWE-KEM.cpp -I/usr/local/include/openfhe -I/usr/local/include/openfhe/pke -I/usr/local/include/openfhe/core -I/usr/local/include/openfhe/cereal -I/usr/local/include/openfhe/binfhe/ -lOPENFHEbinfhe -lOPENFHEcore -lOPENFHEpke -lcrypto -L/usr/local/lib ; ./LWE-KEM

LWE-KEM-compila:
	cd /app/target ; rm LWE-KEM ; g++ -fopenmp -o LWE-KEM LWE-KEM.cpp -I/usr/local/include/openfhe -I/usr/local/include/openfhe/pke -I/usr/local/include/openfhe/core -I/usr/local/include/openfhe/cereal -I/usr/local/include/openfhe/binfhe/ -lOPENFHEbinfhe -lOPENFHEcore -lOPENFHEpke -lcrypto -L/usr/local/lib

# cd /app/target ; rm LWE-KEM ; g++ -o LWE-KEM LWE-KEM.cpp -I/usr/local/include/openfhe -I/# usr/local/include/openfhe/pke -I/usr/local/include/openfhe/core -I/usr/local/include/openfhe/cereal -I/usr/# local/include/openfhe/binfhe/ -lOPENFHEbinfhe -lOPENFHEcore -lOPENFHEpke -lcrypto -L/usr/local/lib

LWE-KEM-loop:
	cd /app/target ; ./loopBash.sh