# Motor Unit Interface (MATLAB)

A MATLAB GUI for post-processing **high-density EMG (HD-EMG) motor unit decomposition** results, built in the **NeuroMechatronics Lab at Carnegie Mellon University** (Prof. Douglas J. Weber).

It takes the output of motor unit decomposition (e.g. from DEMUSE) and lets you inspect and quantify each motor unit.

## What it does

- **Load** a decomposition output (`.mat`) — EMG signals, motor unit pulse trains and force/reference signals (4 kHz sampling)
- **Spatial filtering** — single- or double-differential montage
- **Motor unit action potentials (MUAPs)** — spike-triggered averaging, peak-to-peak and RMS amplitude per channel, plots of single and all MUAPs
- **Recruitment / de-recruitment thresholds** — force level at which each unit turns on and off
- **Conduction velocity** — manual and automatic estimation
- **Spike-train analysis** of the decomposed units

## Repository contents

| File | What it is |
|---|---|
| `MUInterface.m` | Main GUI (`mainGUI`) |
| `Motor_GUI_NML.mlapp` | App Designer version of the GUI |
| `Main_MUTA_ConvexHull.m` | Script version: load data, collect MUAPs, analysis |
| `extractMUAPs_Luigi.m`, `plotAllMUAPs.m`, `plotSingleMUAP.m` | MUAP extraction and plotting helpers |
| `hpbutter.m`, `lpbutter.m` | High- and low-pass Butterworth filters |
| `DEMUSE.zip` | DEMUSE decomposition toolbox used to produce the input files |

## Usage

Open MATLAB, add the folder to the path, update `dataPath` in `MUInterface.m` to your decomposition outputs, and run `MUInterface`.
