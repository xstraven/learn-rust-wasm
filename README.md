# Lorenz Attractor - Rust + WebAssembly Demo

A real-time simulation of the Lorenz attractor implemented in Rust and compiled to WebAssembly for high-performance mathematical computation in the browser.

## Features

- **High-performance computation**: Rust compiled to WebAssembly for near-native performance
- **Real-time visualization**: Smooth 60fps rendering with customizable trail effects
- **Interactive controls**: Adjust Lorenz parameters (σ, ρ, β) and initial conditions in real-time
- **Multiple view modes**: XY, XZ, YZ planes and rotating 3D view
- **Dynamic visualization**: Color-coded trails showing the evolution of the strange attractor

## The Lorenz System

The Lorenz attractor is defined by the system of differential equations:
- dx/dt = σ(y - x)
- dy/dt = x(ρ - z) - y  
- dz/dt = xy - βz

Where σ, ρ, and β are parameters that control the behavior of the system.

## Building and Running

### Prerequisites
- Rust toolchain (rustup)
- wasm-pack (installed automatically by build script)

### Build
```bash
./build.sh
```

### Serve
```bash
python3 -m http.server 8000
```

Then open `http://localhost:8000` in your browser.

## Why Rust + WebAssembly?

This simulation demonstrates the power of Rust compiled to WebAssembly for computational tasks:

1. **Performance**: The simulation performs thousands of floating-point calculations per frame. Rust's zero-cost abstractions provide performance close to native code.

2. **Memory Safety**: Rust's ownership system prevents common bugs like buffer overflows and memory leaks, critical for long-running simulations.

3. **Predictable Performance**: No garbage collection pauses that could cause frame drops during animation.

4. **WebAssembly Benefits**: Runs at near-native speed in the browser while maintaining security through WebAssembly's sandboxed execution environment.

## Controls

- **Start/Stop**: Begin or pause the simulation
- **Reset**: Return to initial conditions and clear the trail
- **Clear Trail**: Remove the trail while keeping the current position
- **Parameters**: Adjust σ (sigma), ρ (rho), β (beta) in real-time
- **Speed**: Control simulation speed (steps per frame)
- **View Modes**: Switch between 2D projections and 3D view
- **Initial Conditions**: Set starting position (X₀, Y₀, Z₀)

## Technical Implementation

- **Rust/WebAssembly**: Core simulation engine with optimized numerical integration
- **Canvas 2D API**: Efficient rendering with alpha blending for trail effects
- **Real-time Parameter Updates**: Modify system behavior during simulation
- **Memory Management**: Configurable point history with automatic cleanup