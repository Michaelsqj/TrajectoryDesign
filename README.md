# TrajectoryDesign

`Simultaneous angiography & perfusion of ASL`

### Base Cone Types

- Radial
- Johnson
- Gurney
- other variants

### Metrics:

- Density - radius
- SNR efficiency
- Histogram of weights distribution
- FWHM


Entry points are two functions `main_angi.m` & `main_perf.m`.

### Process

1. Generate ktraj
    1. Set all parameters
    2. Generate trajectory
        1. Rotate
            1. Generate base trajectory (`base_g`, `base_k`)
            2. Rotate base trajectory
        2. Stack
            1. Generate base trajectory (`gen_StackCones.m`)
            2. Scale & rotate base trajectory
    4. Calculate trajectory intrinsic metrics
2. Reconstruct phantom images
3. Save and print

***TODO***
    - Test `main_angi.m`
    - Implement `main_perf.m`
    - Implement stack of cones
    - Implement full spoke symmetric cone