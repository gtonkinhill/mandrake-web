CXX = em++
# SOURCES = pathosce_bindings.cpp

SOURCES= pp-src/dist/dist.o pp-src/dist/matrix_ops.o pp-src/reference.o pp-src/sketch/seqio.o pp-src/sketch/sketch.o pp-src/sketch/countmin.o pp-src/dist/linear_regression.o pp-src/random/rng.o pp-src/random/random_match.o pathosce_bindings.cpp


# CXXFLAGS+=-Wall -Wextra -std=c++14 -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -fPIC 


# CXXFLAGS += -O3 -g\
# 	--closure 1 \
# 	-flto -s STRICT=1 \
# 	-s MODULARIZE=1 \
# 	-s WASM=1 \
# 	-s NO_DISABLE_EXCEPTION_CATCHING \
# 	-s USE_BOOST_HEADERS=1 \
# 	-s USE_ZLIB=1 \
# 	-s FORCE_FILESYSTEM=1 \
# 	-s ALLOW_MEMORY_GROWTH=1 \
# 	-lworkerfs.js \
# 	--bind \
# 	# -DNOEXCEPT \
# 	# -DJSON_NOEXCEPTION \
# 	# -fno-exceptions \
# 	# -s LLD_REPORT_UNDEFINED=1 \
# 	# -s ERROR_ON_UNDEFINED_SYMBOLS=0

CXXFLAGS = -O3 -g \
	-s WASM=1 \
	-s NO_DISABLE_EXCEPTION_CATCHING \
	-s USE_BOOST_HEADERS=1 \
	-s USE_ZLIB=1 \
	-s FORCE_FILESYSTEM=1 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-lworkerfs.js \
	--bind \
	-s LLD_REPORT_UNDEFINED \
	--no-entry

# CXXFLAGS = -O3 -g  \
# 		-s WASM=1 \
# 				-flto --bind -s STRICT=1 \
# 				-s ALLOW_MEMORY_GROWTH=1 \
# 				-s USE_ZLIB=1 \
# 				-s USE_BOOST_HEADERS=1 \
# 				-s FORCE_FILESYSTEM=1 \
# 				-s MODULARIZE=1 \
# 				-lworkerfs.js \
# 				--bind \
# 				# --closure 1 \
# 				# -fopenmp -lomp \
# 				# -s "EXPORTED_FUNCTIONS=['_malloc']" \
# 				# -s 'EXTRA_EXPORTED_RUNTIME_METHODS=["FS"]' \
# 				# -s EXPORT_NAME=WebSketch \
# 				-Wall -Wextra -std=c++14 \
# 				# -fno-exceptions \
# 				# -DNOEXCEPT \
# 				# -DJSON_NOEXCEPTION \
# 				# -s DISABLE_EXCEPTION_CATCHING=1 \

# CXXFLAGS = -O3 -s ASSERTIONS=1  \
# 				--closure 1 \
# 				-flto --bind -s STRICT=1 \
# 				-s ALLOW_MEMORY_GROWTH=1 \
# 				-s USE_ZLIB=1 \
# 				-s USE_BOOST_HEADERS=1 \
# 				-s FORCE_FILESYSTEM=1 \
# 				-s MODULARIZE=1 \
# 				-s "EXPORTED_FUNCTIONS=['_malloc']" \
# 				-s 'EXTRA_EXPORTED_RUNTIME_METHODS=["FS"]' \
# 				-lworkerfs.js \
# 				-fopenmp -lomp \
# 				# -s EXPORT_NAME=WebSketch \
# 				-Wall -Wextra -std=c++14 \
# 				# -fno-exceptions \
# 				# -DNOEXCEPT \
# 				# -DJSON_NOEXCEPTION \
# 				# -s DISABLE_EXCEPTION_CATCHING=1 \

# LDLIBS = -lz -lhdf5_cpp -lhdf5 -lopenblas -llapack -lgomp -lomp -lgfortran
# LDLIBS += -lpthread -lgfortran -lm -ldl -lrt -lomp
# LDFLAGS=-Wl,-as-needed -lomp -pthread -lomp
# -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include

LIBLOC = ${CONDA_PREFIX}
CXXFLAGS += -I"." -I'./pp-src' -I'./armadillo-10.6.2/include' -I$(LIBLOC)/include/eigen3 -I$(LIBLOC)/include/highfive -I'./armadillo-10.6.2' -I/usr/local/opt/libomp/include 
# CXXFLAGS += -DARMA_DONT_USE_WRAPPER
LDFLAGS = -L/usr/local/opt/libomp/lib -L'./armadillo-10.6.2/lib'
# LDLIBS = -larmadillo -lz -lhdf5_cpp -lhdf5 -lopenblas -llapack -lgomp -lomp -lgfortran

# LDLIBS = -lz -lhdf5_cpp -lhdf5 -lopenblas -llapack -lgomp -lomp -lgfortran -lpthread -lgfortran -lm -ldl -lrt -lomp
# LDFLAGS=-Wl,-as-needed -lomp -L'./armadillo-10.6.2/'

# web: CXX = em++
# # optimised compile options
# # NB turn exceptions back on for testing
# web: CXXFLAGS = -O3 -s ASSERTIONS=1  \
# 				--closure 1 \
# 				-flto --bind -s STRICT=1 \
# 				-s ALLOW_MEMORY_GROWTH=1 \
# 				-s USE_ZLIB=1 \
# 				-s MODULARIZE=1 \
# 				-s "EXPORTED_FUNCTIONS=['_malloc']" \
# 				-s 'EXTRA_EXPORTED_RUNTIME_METHODS=["FS"]' \
# 				-s EXPORT_NAME=WebSketch \
# 				-Wall -Wextra -std=c++14 \
# 				# -fno-exceptions \
# 				# -DNOEXCEPT \
# 				# -DJSON_NOEXCEPTION \
# 				# -s DISABLE_EXCEPTION_CATCHING=1 \
# web: CPPFLAGS += -DWEB_SKETCH
# web: LDFLAGS = -lworkerfs.js

# web: ./pp-src/dist/dist.o ./pp-src/dist/matrix_ops.o ./pp-src/reference.o ./pp-src/sketch/seqio.o ./pp-src/sketch/sketch.o ./pp-src/database/database.o ./pp-src/sketch/countmin.o ./pp-src/dist/linear_regression.o ./pp-src/random/rng.o ./pp-src/random/random_match.o pathosce_bindings.cpp
# 	$(LINK.cpp) $^ -o pathosce.js 


build: $(SOURCES)
	$(LINK.cpp)  $^ -o ../static/js/pathosce.js

# build: $(SOURCES)
# 	$(EMCC) $(EMFLAGS) $(SOURCES) $(LDFLAGS) $(CPPFLAGS)  -o ../static/js/pathosce.js  $(LDLIBS)

clean :
	rm ../static/js/pathosce.js ../static/js/pathosce.wasm ./*.o ./pp-src/*.o ./pp-src/*/*.o