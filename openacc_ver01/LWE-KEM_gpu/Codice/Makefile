CXX = nvc++ 
# CXXFLAGS = -std=c++20 -O2 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/include/botan-3 -I/opt/homebrew/opt/openssl@3/include
# LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp -L/opt/homebrew/lib -lbotan-3 -L/opt/homebrew/opt/openssl@3/lib -lcrypto
CXXFLAGS = -std=c++20 -O2 -acc -gpu=cc90,managed -Minfo=accel
LDFLAGS  = -lcrypto #-fopenmp


all: lwe_kem

clean: 
	rm -f lwe_kem

lwe_kem: Codice_modificato.cpp
#	$(CXX) $(CXXFLAGS) Codice_modificato.cpp pke.cpp utils.cpp hash.cpp hash_openssl.cpp noise.cpp kem.cpp $(LDFLAGS) -o lwe_kem
	$(CXX) $(CXXFLAGS) Codice_modificato.cpp pke.cpp utils.cpp hash_openssl.cpp noise.cpp kem.cpp $(LDFLAGS) -o lwe_kem

# NEW PART
# --- bench params ---
#KEYS ?= 50
#OUT  ?= bench.csv

#N_LIST = 512 1024 2048 4096 8192
#T_LIST = $(shell seq 2 2 36)
#T_LIST = 1 2 4 6 8 12 16 24 32 36

#run: lwe_kem
#	OMP_DYNAMIC=FALSE OMP_PROC_BIND=true OMP_PLACES=cores OMP_NUM_THREADS=$(THREADS) ./lwe_kem

#bench_arm: lwe_kem
#	@echo "n;threads;avg_keygen_us;avg_encaps_us;avg_decaps_us;total_s;mismatches" > $(OUT)
#	@for n in $(N_LIST); do \
	  for t in $(T_LIST); do \
	    echo "Running n=$$n threads=$$t"; \
	    OMP_DYNAMIC=FALSE OMP_PROC_BIND=true OMP_PLACES=cores OMP_NUM_THREADS=$$t ./lwe_kem $$n $(KEYS) >> $(OUT); \
	  done; \
	done