#!/bin/bash

echo "Building Rust WebAssembly module..."

if ! command -v wasm-pack &> /dev/null; then
    echo "wasm-pack not found. Installing..."
    curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
fi

wasm-pack build --target web --out-dir pkg --no-typescript

if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "You can now serve the project with: python3 -m http.server 8000"
    echo "Then open http://localhost:8000 in your browser"
else
    echo "Build failed!"
    exit 1
fi