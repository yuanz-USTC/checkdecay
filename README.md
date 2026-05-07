# checkdecay

ROOT macro project for estimating the J/psi -> e+e- PID efficiency vs J/psi `pT` and rapidity.

## Files

- `JpsiDecay.C`: standalone J/psi decay utilities. It can generate uniform J/psi `pT/y/phi` and decay isotropically to dielectrons.
- `PidEfficiencyMaps.C`: reads the three leg PID efficiency maps vs `p` and `eta`.
- `run_checkdecay.C`: main macro. It supports generated decays and future stored tree input.
- `input/`: put PID maps and optional J/psi -> ee tree here.
- `output/`: output ROOT files are written here.

## PID Map Convention

The macro expects three `TH2` maps with x-axis = leg momentum `p` and y-axis = leg `eta`:

- `hPidEffAvg`: average leg efficiency, independent of whether the other leg passes.
- `hPidEffIfOtherPass`: remaining-leg efficiency if the other leg passes PID.
- `hPidEffIfOtherFail`: remaining-leg efficiency if the other leg fails PID.

The names are configurable in `run_checkdecay`.

## Run

From this directory:

```bash
  root -l -q 'run_checkdecay.C("simulate",1000000)'
```

Default output:

```text
output/jpsi_pid_efficiency.root
```

The main result histogram is:

```text
hJpsiPidEfficiency
```

## Tree Mode

When the input tree is ready, put it under `input/` and run:

```bash
  root -l -q 'run_checkdecay.C("tree",1000000,"input/pid_efficiency_maps.root","output/jpsi_pid_efficiency.root",80,0,20,80,-4,4,"input/jpsi_to_ee_tree.root","jpsiTree")'
```

Required branch/leaf names:

- `jpsi_pt`, `jpsi_eta`, `jpsi_phi`, `jpsi_energy`
- `leg1_pt`, `leg1_eta`, `leg1_phi`
- `leg2_pt`, `leg2_eta`, `leg2_phi`

## Defaults Before Inputs Exist

If `input/pid_efficiency_maps.root` is not present, the macro still runs with fallback efficiencies:

- average: `0.85`
- if other leg passes: `0.90`
- if other leg fails: `0.70`

These are only placeholders so the workflow can be tested before the real maps are available.
