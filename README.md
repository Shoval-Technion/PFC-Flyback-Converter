# Arduino-Oriented PFC Flyback Converter

**Characterization, Analysis and Verification of a PFC Flyback Converter**  
PEARL Power Electronics Laboratory | Spring 2026

Project Team: Anat Rom, Shoval Buzaglo

Academic Supervisor: Eric Herbelin


## Project Overview

This academic project investigates a **12 V / 9 W isolated Flyback converter** with power factor correction (PFC). The goal is to evaluate whether a control strategy suitable for an Arduino-class microcontroller can replace a dedicated analog PFC controller while maintaining high input power factor under different output-load conditions.

The work began with a study of the **EV44010 / MP44010** reference converter. Two LTspice models of our converter were then developed and compared:

- **Without PFC:** a conventional fixed-frequency, voltage-mode PI controller with a large rectified-input capacitor.
- **With PFC:** fixed-frequency, discontinuous-conduction-mode (DCM) operation with a slowly varying constant-on-time command and a small input capacitor.

The report documents **simulation results**. An Arduino-controlled hardware prototype and laboratory validation are future steps; the reported results should not be interpreted as hardware measurements.

## Design Specifications

| Parameter | Design value |
| --- | --- |
| Input | 35.36 V RMS AC, 50 Hz (50 V peak) |
| Output | 12 V DC, nominal 9 W |
| Converter | Isolated, single-switch Flyback |
| Switching frequency | 100 kHz, fixed |
| Conduction mode | DCM |
| Primary inductance | 90 µH |
| Secondary inductance | 12.7 µH |
| PFC control | Slowly varying constant on-time |
| Nominal on-time | ~3.43 µs |
| Feedback-filter corner | ~10.6 Hz |

## Control Concept

For DCM operation, the primary current starts near zero at the beginning of each switching cycle. With a nearly constant switch on-time, the cycle-averaged input current is approximately

```text
I_in,avg(t) = V_in(t) * T_on^2 / (2 * L_p * T_s)
```

Thus, with nearly fixed `T_on` and switching period `T_s`, the average input current is approximately proportional to the instantaneous rectified input voltage. The converter behaves approximately like an emulated input resistance. A deliberately slow voltage-feedback path adjusts the **average power demand** without strongly following the 100 Hz output ripple.

Key changes between the two models include reducing the rectified-input capacitor from **125 µF to 10 nF** and changing the voltage-regulation dynamics to maintain an almost constant on-time during a line cycle.

## Simulation and Analysis

- **LTspice:** circuit models, transient simulation, switching waveforms, and the behavioral representation of the control algorithm.
- **MATLAB:** post-processing the simulated voltage and current waveforms, including calculations of power factor (PF), total harmonic distortion (THD), power, efficiency, and comparative plots.
- **Analysis interval:** 0.46–0.50 s for the reported steady-state comparisons.
- **Load cases:** resistive, inductive, capacitive, and repeated load steps.

**Analysis note:** The reported analysis was performed using MATLAB scripts, rather than treating LTspice `.meas` directives as the source of the presented calculations. Reproduction of the results requires the corresponding MATLAB scripts and waveform data/export procedure.

## Selected Simulation Results

| Test case | Power factor | Input-current THD | Efficiency |
| --- | ---: | ---: | ---: |
| Without PFC (resistive baseline) | 0.5698 | 139.77% | 80.13% |
| PFC — resistive load | 0.9976 | 6.95% | 85.12% |
| PFC — inductive load | 0.9975 | 7.01% | 85.07% |
| PFC — capacitive load | 0.9975 | 6.97% | 85.59% |

The PFC model maintained a power factor above 0.99 for the tested load-step cases up to 100% → 200%; the extreme 100% → 400% step reduced PF to approximately **0.8846**. Results depend on the modeled operating points and are not yet hardware-validated.

## Repository Files

The repository is intended to contain:

```text
LTspice/         Circuit schematics, netlists and required device models
MATLAB/          Analysis scripts and instructions
Documentation/   Final report and presentation
References/      Technical documentation and supporting literature
```

## Future Work

- Implement the slow control loop and timer-based 100 kHz PWM on an Arduino-class microcontroller.
- Implement and characterize isolated voltage feedback and fast hardware current limiting.
- Build a prototype, verify DCM operation, and compare laboratory measurements against simulations.
- Investigate output-voltage transients and PFC performance under larger dynamic load steps.

## Documentation

See the project final report, **“Characterization, Analysis and Verification of a PFC Flyback Converter”**, for the theoretical background, circuit models, methodology, and full comparison of simulation cases.
