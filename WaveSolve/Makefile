clean:
	rm -rf out

build-tests:
	mkdir -p out
	mkdir -p out/test
	mkdir -p out/test/wasm
	emcc -fno-exceptions -c tests/main.cpp -o out/test/wasm/main.wasm
	emcc -fno-exceptions -c tests/fft.test.cpp -o out/test/wasm/fft.test.wasm
	emcc -fno-exceptions -c src/fft.cpp -o out/test/wasm/fft.wasm
	emcc -fno-exceptions out/test/wasm/main.wasm out/test/wasm/fft.test.wasm out/test/wasm/fft.wasm -o out/test/wasm/test.wasm

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
