# single_leg_efficiency

Extract leg PID efficiency maps from reconstructed J/psi candidates.

Input default:

```text
../input/LHC25b16_mcrecgen.root
DF_2363808317891232/O2rtdiemtreerec
```

Event preselection:

```text
(fFilterMap & (1u << 0)) != 0
fIsCorrectAssoc == true
```

Leg PID cut:

```text
-2.5 <= fTPCNSigmaEl <= 3.0
fTPCNSigmaPi > 3.0
fTPCNSigmaPr > 3.0
```

Run from this directory:

```bash
  root -l -b -q 'make_single_leg_efficiency_maps.C'
```

Default output:

```text
../input/pid_efficiency_maps.root
```

The output names match `../PidEfficiencyMaps.C`:

- `hPidEffAvg`
- `hPidEffIfOtherPass`
- `hPidEffIfOtherFail`

The file also contains numerator and denominator histograms for each map.
