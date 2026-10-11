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
root -l -q 'run_checkdecay.C("tree",1000000,"input/pid_efficiency_maps.root","output/jpsi_pid_efficiency.root",80,0,20,80,-4,4,"input/jpsi_to_ee_tree.root","jpsiTree",12345,0.85,0.90,0.70,"hPidEffAvg","hPidEffIfOtherPass","hPidEffIfOtherFail","gen",0)'
```

The last arguments select the tree level and cuts:

- `treeLevel = "gen"` or `"rec"`
- `mcDecisionBit = 0` by default; set to `-1` to disable this cut.
- `filterMapBit = 0` by default for rec trees; set to `-1` to disable this cut.
- `requireCorrectAssoc = true` by default for rec trees.

### Gen-Level Tree

Required branch/leaf names:

- J/psi: `fE`, `fPt`, `fEta`, `fPhi`
- leg 1: `fPtMC1`, `fEtaMC1`, `fPhiMC1`
- leg 2: `fPtMC2`, `fEtaMC2`, `fPhiMC2`
- selection: `fMcDecision`

Default selection:

```text
(fMcDecision & (1u << 0)) != 0
```

### Rec-Level Tree

Example:

```bash
root -l -q 'run_checkdecay.C("tree",1000000,"input/pid_efficiency_maps.root","output/jpsi_pid_efficiency_rec.root",80,0,20,80,-4,4,"input/jpsi_to_ee_tree.root","jpsiTree",12345,0.85,0.90,0.70,"hPidEffAvg","hPidEffIfOtherPass","hPidEffIfOtherFail","rec",0,0,true)'
```

Required branch/leaf names:

- J/psi: `fMass`, `fPt`, `fEta`, `fPhi`
- leg 1: `fPt1`, `fEta1`, `fPhi1`
- leg 2: `fPt2`, `fEta2`, `fPhi2`
- selection: `fMcDecision`, `fFilterMap`, `fIsCorrectAssoc`

Default selection:

```text
(fMcDecision & (1u << 0)) != 0
(fFilterMap & (1u << 0)) != 0
fIsCorrectAssoc == true
```

## Defaults Before Inputs Exist

If `input/pid_efficiency_maps.root` is not present, the macro still runs with fallback efficiencies:

- average: `0.85`
- if other leg passes: `0.90`
- if other leg fails: `0.70`

 These are only placeholders so the workflow can be tested before the real maps are available.
