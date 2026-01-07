CXX = g++
CXXFLAGS = -std=c++20 -O2 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/include/botan-3 -I/opt/homebrew/opt/openssl@3/include
LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp -L/opt/homebrew/lib -lbotan-3 -L/opt/homebrew/opt/openssl@3/lib -lcrypto

all: lwe_kem

clean: 
	rm -f lwe_kem

lwe_kem: Codice_modificato.cpp
	$(CXX) $(CXXFLAGS) Codice_modificato.cpp pke.cpp utils.cpp hash.cpp hash_openssl.cpp noise.cpp kem.cpp $(LDFLAGS) -o lwe_kem