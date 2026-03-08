CXX      = nvc++
CXXFLAGS = -std=c++20 -O2 -acc -gpu=cc80,managed -Minfo=accel
LDFLAGS  = -lcrypto

SRCS = Codice_modificato.cpp pke.cpp utils.cpp hash_openssl.cpp noise.cpp kem.cpp

all: lwe_kem

clean:
	rm -f lwe_kem

lwe_kem: $(SRCS) xorshift.h
	$(CXX) $(CXXFLAGS) $(SRCS) $(LDFLAGS) -o lwe_kem

# ---------- benchmark ----------
N_LIST = 10
DIM_LIST = 512 1024 2048 4096

bench: lwe_kem
	@echo "n;avg_keygen_us;avg_encaps_us;avg_decaps_us;total_s;mismatches"
	@for d in $(DIM_LIST); do \
	  echo ">>> n=$$d"; \
	  ./lwe_kem $(N_LIST) $$d; \
	done
