# unstable-invasion
Code for simulating invasion of unstable populations

This code simulates the invasion of unstable populations into stable "healthy" populations.
The simulations are performed on 2D & 3D flat and radial geometries.

## Dependencies
- Boost
- MPI with C++ bindings

### macOS
```
$ brew install boost
$ brew edit open-mpi
```
Add the following line to the args list, beginning on line 55:
```
--enable-mpi-cxx
```

## Building

```
mkdir -p build
cd build
cmake ..
make
```

Compiled binaries can then be found in `build/bin`
