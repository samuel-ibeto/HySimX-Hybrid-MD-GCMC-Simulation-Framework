# HySimX: Hybrid MD/GCMC Simulation Framework
This is a Hybrid Molecular Dynamics-Grand Canonical Monte Carlo Simulation Framework for eXtended Adsorption and Structural Relaxation in Host-Guest Systems. A Python-based automation framework that orchestrates iterative cycles of gas insertion and structural relaxation within a single, reproducible workflow using GROMACS and CASSANDRA softwares.[link to article](PUT URL HERE WHEN PUBLISHED)


## **Dependencies**

### **1. Python**
- **Python 3.10.12** (https://www.python.org/downloads/release/python-31012/)
  
### **2. GROMACS**
- Tested with **GROMACS 2022.4** or higher (https://www.gromacs.org/)

### **3. Cassandra**
- **Cassandra v 1.3.1**  (https://github.com/MaginnGroup/Cassandra)

### **4. Open Babel**
- **Version >= 3.1** (https://github.com/openbabel/openbabel)
---
## **Overview**

HySimX couples **molecular dynamics (MD)** simulations in **GROMACS** with **grand canonical Monte Carlo (GCMC)** sampling in **Cassandra** within a single reproducible workflow.

The framework executes iterative **adsorption–relaxation cycles**, combining:

- **MC-based molecule insertion/deletion**
- **MD-based structural equilibration**

This enables the realistic modeling of:

- **adsorption**
- **transport**
- **structural adaptation**
- **host–guest interactions**
- **free-volume evolution**

in flexible materials such as **polymers**, **porous frameworks**, **liquids**, and **heterogeneous interfaces**.

## **Key Features**

- **Automated hybrid MD-GCMC workflow**
- **Iterative adsorption-relaxation cycles**
- **Compatible with GROMACS and Cassandra**
- **Modular and extensible design**
- **Topology and coordinate updates between cycles**
- **Restart and checkpoint functionality**
- **Suitable for high-throughput simulations**
- **Designed for diverse host-guest systems**
---

## **Directory Structure**

```text
HySimX/
├── input/                          # Initial configuration files
│   ├── start_conf.gro
│   └── start_conf.tpr
│
├── cassandra/                      # Cassandra-related files
│   ├── speciesX/                   # Species folders
│   └── cass_inp.inp                # Cassandra input file(s)
│   └── *.mcf                       # Cassandra molecular configuration file(s)
│
├── gromacs/                        # GROMACS-related files
│   ├── system.top                  # Combined topology file
│   ├── min.mdp                     # Energy minimization parameters
│   └── md.mdp                      # MD relaxation parameters (NVT or NPT)
│
├── scripts/                        # Utility scripts
│   └── xyz_to_gro.py               # Converts Cassandra .xyz -> GROMACS .gro
│
├── cycle_XXX/                      # Auto-generated cycle directories
│   ├── cass_inp.inp
│   ├── output-XXX.gro
│   ├── *.prp
│   └── xyz_to_gro.py
│
├── config.yml                      # Main configuration file
└── main.py                         # Main workflow driver
```
---

## **Configuration File (config.yml)**

The `config.yml` file defines the key parameters required to execute the HySimX workflow. These include:

- **Working directory:** `work_dir`
- **Number of simulation cycles:** `cycles`
- **Initial structure files:** starting `.gro` and `.tpr` files
- **Executable paths:** `gmx`, `cassandra`, `obabel`, and `python`
- **Cassandra settings:** input file name, species folder prefix, completion string, and `.prp` output pattern
- **GROMACS settings:** topology file, minimization `.mdp`, and MD relaxation `.mdp` files
- **Molecule mapping:** identifiers for adsorbents and gas species
- **Backup options**
---


## **Usage**

This section provides a brief overview of how to use the HySimX framework. A working knowledge of GROMACS and Cassandra  is assumed. Users are encouraged to consult the official documentation for both packages for detailed guidance on input file preparation, simulation settings, and parameter selection.

Before launching extended production runs, it is strongly recommended to perform a short single-cycle test run to verify that the input parameters, file paths, and workflow settings are correctly configured. This is particularly important for complex or multicomponent systems where small setup inconsistencies may lead to errors or unintended simulation behavior.

### **1. Edit the configuration file**
Update the simulation settings in `config.yml` before starting a run.

```bash
nano config.yml
```
Next make sure all the required 
### **2. Run the workflow**
To start a new simulation workflow, run:
```bash
python3 main.py config.yml
```
### **3. Monitor the outputs**
During execution, each cycle generates a cycle_XXX directory containing the corresponding Cassandra and GROMACS outputs.

## **Checkpointing**
After each completed cycle, a `checkpoint.json` file is automatically written to the working directory. This file stores:
- the most recently completed cycle
- elapsed simulation time
- box information
- molecule counts
- paths to intermediate files

The complete cycle history is preserved to support traceability, restart safety, and reproducibility.

## **Automatic Restart Detection**
HySimX includes a robust restart system that enables interrupted simulations to be resumed efficiently and completed workflows to be extended seamlessly. When a restart is requested, the framework automatically identifies the last fully completed cycle in the working directory and resumes from the next cycle. This prevents unnecessary repetition of completed work.

### **Restart Options**
Use the command-line flag `--restart` to resume or extend a previous run:
```bash
python3 main.py config.yml --restart
```
- If no number is provided, the framework will continue until the YAML-specified total cycles are completed.  
- If a number is provided, e.g., `--restart 5`, the framework will run that many additional cycles beyond the last completed cycle.  

### **Restart Example**
**Example 1**

Last fully completed cycle: 7

YAML cycles specified: 10

Command:
```bash
 python3 main.py config.yml --restart
```
Result:

Cycles 8, 9, and 10 will be executed.

**Example 2**

Last fully completed cycle: 8

YAML cycles specified: 10

Command:
```bash
python3 main.py config.yml --restart 5
```

Result: 

Cycles 9, 10, 11, 12 and 13 will be executed.

- In this example, 5 indicates that five additional cycles will be executed beyond the last completed cycle.

