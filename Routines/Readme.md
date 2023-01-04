## Processes

Main routines used by simulations.

- `DefineStructs.jl` (struct)
    <br> Defines the structs created during the simulations: 
    - `ObsDetect.jl` contains the observed epidemic (e.g., observed time series of cases), 
    - `SimInfect.jl` contains the simulated undetected epidemic (e.g., the new infections, the final size of the epidemic and the secondary infections) and 
    - `SimDetect.jl` contains the detected infections in our simulations (e.g., the time series of the simulated cases, the per-day cumulative number of cases).

- `InfectionProcess.jl` (function)
    <bf> Models disease transmission starting from a single individual. Requires the context-specific parametrization and returns a `SimInfect` struct.

- `DetectionProcess.jl` (function)
    <bf> Models the detection of infected individuals. Requires the object `SimInfect` and the final size of the observed data, `ObsDetect`, and returns a `SimDetect` struct.

- `DetectionProcess.jl` (function)
    <bf> Routine to display all results, which are summarized using median values and 95% interquantile ranges (95%IqR; values between the 2.5th and the 97.5th percentiles).

### References
To learn about `struct`s: <a href="https://docs.julialang.org/en/v1/base/base/#struct" rel="_blank">https://docs.julialang.org/en/v1/base/base/#struct</a>
