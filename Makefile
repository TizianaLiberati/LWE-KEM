CXX = g++
CXXFLAGS = -std=c++20 -O2 -I/opt/homebrew/include/botan-3
LDFLAGS = -L/opt/homebrew/lib -lbotan-3

all: lwe_kem

lwe_kem: Codice_modificato.cpp
	$(CXX) $(CXXFLAGS) Codice_modificato.cpp $(LDFLAGS) -o lwe_kem