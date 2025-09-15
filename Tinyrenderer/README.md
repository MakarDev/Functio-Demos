# TinyRenderer Optimizations

This folder contains both the **original** and **optimized** C++ versions of selected code from [ssloy/tinyrenderer](https://github.com/ssloy/tinyrenderer).

For clarity and ease of comparison, each variant has been consolidated into a **single `.cpp` file** rather than being split across multiple headers and sources.  
This makes it straightforward to diff, benchmark, and study the changes.

## Contents
- `original_*.cpp` – Direct copies of TinyRenderer’s original implementations.
- `optimized_*.cpp` – Modified versions with performance-oriented changes (vectorization, cache-friendly loops, etc.).

## Notes
- These files are self-contained and meant for demonstration and performance benchmarking.
- The goal is to showcase optimization techniques, not to replace TinyRenderer’s structure.
- Build instructions remain simple: just compile each `.cpp` independently with your compiler of choice.
