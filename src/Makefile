EMCC = em++
SOURCES = pathosce_bindings.cpp

EMFLAGS = -O3 -g \
	-s WASM=1 \
	-s NO_DISABLE_EXCEPTION_CATCHING \
	-s USE_BOOST_HEADERS=1 \
	-s USE_ZLIB=1 \
	-s FORCE_FILESYSTEM=1 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-lworkerfs.js \
	--bind

build: $(SOURCES)
	$(EMCC) $(EMFLAGS) $(SOURCES) -o ../static/js/pathosce.js

clean :
	rm pathosce.js pathosce.wasm