clean:
	rm -rf out

build-tests:
	mkdir -p out
	mkdir -p out/test
	mkdir -p out/test/wasm
	emcc -fexceptions -s SIDE_MODULE=1 -s STRICT=1 -c tests/main.cpp -o out/test/wasm/main.o
	emcc -fexceptions -s SIDE_MODULE=1 -s STRICT=1 -c tests/fft.test.cpp -o out/test/wasm/fft.test.o
	emcc -fexceptions -s SIDE_MODULE=1 -s STRICT=1 -c src/fft.cpp -o out/test/wasm/fft.o
	emcc -fexceptions -s MAIN_MODULE=1 -s STRICT=1 -s STANDALONE_WASM=1 out/test/wasm/main.o out/test/wasm/fft.test.o out/test/wasm/fft.o -o out/test/wasm/test.wasm

test:
	wasmtime out/test/wasm/test.wasm

build:
	mkdir -p out
	mkdir -p out/build
	mkdir -p out/build/wasm

dist:
	mkdir -p out
	mkdir -p out/dist
	mkdir -p out/dist/wasm
	mkdir -p out/dist/asmjs
