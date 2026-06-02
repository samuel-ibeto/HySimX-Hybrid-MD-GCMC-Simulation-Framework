#!/usr/bin/env python3
"""
main.py

Automates GCMC (Cassandra) <-> MD (GROMACS) cycles 

Run:
    python main.py config.yml
"""

import sys
import subprocess
import shutil
import time
import json
from pathlib import Path
from datetime import datetime
import yaml
import re
import os

# ------------------------------
# Utility helpers
# ------------------------------
def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def timenow():
    return time.perf_counter()

def run_cmd(cmd, cwd=None, capture_output=False, check=True, shell=False):
    """Run command using subprocess.run; raise RuntimeError on failure (if check True)"""
    print(f"[{now()}] RUN: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    try:
        if capture_output:
            proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=shell)
            out = proc.stdout
            if check and proc.returncode != 0:
                raise RuntimeError(f"Command failed (rc {proc.returncode})\nOutput:\n{out}")
            return proc.returncode, out
        else:
            proc = subprocess.run(cmd, cwd=cwd, shell=shell)
            if check and proc.returncode != 0:
                raise RuntimeError(f"Command failed (rc {proc.returncode})")
            return proc.returncode, None
    except FileNotFoundError as e:
        raise RuntimeError(f"Executable not found: {e}")

def backup_file(path: Path):
    if not path.exists():
        return None
    ts = datetime.now().strftime("%Y%m%d%H%M%S")
    bak = path.with_suffix(path.suffix + f".bak_{ts}")
    shutil.copy2(path, bak)
    print(f"[{now()}] Backed up {path} -> {bak}")
    return bak

# ------------------------------
# Config & validation
# ------------------------------
def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def ensure_executable(path_str):
    p = Path(path_str)
    if p.exists() and os.access(str(p), os.X_OK):
        return True
    # try which
    from shutil import which
    if which(path_str):
        return True
    return False

def validate_start_files(cfg, work_dir: Path):
    print(f"[{now()}] Validating required files and folders...")
    # required files list
    required = []
    required.append(work_dir / cfg['start']['gro'])
    required.append(work_dir / cfg['start']['tpr'])
    required.append(work_dir / cfg['cassandra']['cass_inp'])
    required.append(work_dir / cfg['gromacs']['topology'])
    required.append(work_dir / cfg['gromacs']['min_mdp'])
    required.append(work_dir / cfg['gromacs']['md_mdp'])
    required.append(work_dir / cfg['xyz_to_gro']['script'])
    # check mcf files listed in ALL.inp
    cass_inp = work_dir / cfg['cassandra']['cass_inp']
    if not cass_inp.exists():
        raise RuntimeError(f"Cassandra input file not found: {cass_inp}")
    mol_files = parse_molecule_files_from_cass_inp(cass_inp)
    
    # verify species folder(s)  are present 
    species_prefix = cfg['cassandra'].get('species_folder_prefix', 'species')
    
    species_dirs = sorted([
        p for p in cass_inp.parent.iterdir()
        if p.is_dir() and p.name.startswith(species_prefix)
    ])   
    if not species_dirs:   
        raise RuntimeError( 
            f"No species folders found in {cass_inp.parent}. "
            f"Expected directories starting with '{species_prefix}'."         
        )                
    print(f"[{now()}] Found species folders: {[d.name for d in species_dirs]}")   
    
    #here the cassandra directory is defined 
    cass_dir = (work_dir / cfg['cassandra']['cass_inp']).parent
    
    # check presence and collect missing mcf
    missing = []
    for fname, cnt in mol_files:
        # if mapping provided, ensure mapping matches  
        # first check in cassandra folder then work_dir
        cand1 = cass_dir / fname
        cand2 = work_dir / fname
        if not cand1.exists() and not cand2.exists():
            missing.append(fname)
    if missing:
        raise RuntimeError(f"Missing .mcf files required by Cassandra input: {missing}")
    # validate executables
    exe_cfg = cfg['executables']
    for key in ['gmx', 'cassandra', 'obabel', 'python']:
        if not ensure_executable(exe_cfg.get(key)):
            print(f"[{now()}] WARNING: Executable for '{key}' not found or not executable: {exe_cfg.get(key)}. Ensure path is correct.")
    print(f"[{now()}] Validation ok.")
    return mol_files

# ------------------------------
# File parsers & editors
# ------------------------------
def gro_box_size_nm(gro_path: Path):
    """Read .gro last line and parse numeric box dims; return average (nm) scalar"""
    txt = gro_path.read_text().strip().splitlines()
    if not txt:
        raise RuntimeError(f"{gro_path} is empty")
    last = txt[-1].strip()
    parts = last.split()
    nums = []
    for p in parts:
        try:
            nums.append(float(p))
        except:
            pass
    if not nums:
        raise RuntimeError(f"No numeric box found in {gro_path}")
    if len(nums) >= 3:
        avg = sum(nums[:3]) / 3.0
        return avg
    return nums[0]

def edit_cassandra_inp_replace_box_and_start(cass_inp_path: Path, box_angstrom: float, start_xyz_name: str, backup=False):
    """Modify Cassandra input file in place: replace # Box_Info cubic value and # Start_Type filename (.xyz)"""
    text = cass_inp_path.read_text()
    lines = text.splitlines()
    # find # Box_Info
    bi = None
    for i,ln in enumerate(lines):
        if ln.strip().startswith("# Box_Info"):
            bi = i
            break
    if bi is None:
        raise RuntimeError("Cannot find '# Box_Info' in the Cassanda input file")
    if bi + 3 >= len(lines):
        raise RuntimeError("Unexpected Cassanda input format near # Box_Info")
    old_val = lines[bi+3]
    lines[bi+3] = f"{float(box_angstrom):.6f}"
    # find # Start_Type
    si = None
    for i,ln in enumerate(lines):
        if ln.strip().startswith("# Start_Type"):
            si = i
            break
    if si is None:
        raise RuntimeError("Cannot find '# Start_Type' in Cassanda input")
    # find next non-empty config line
    cfg_idx = None
    for j in range(si+1, si+6):
        if j < len(lines) and lines[j].strip() != "":
            cfg_idx = j
            break
    if cfg_idx is None:
        raise RuntimeError("Couldn't find Start_Type config line")
    cfg_line = lines[cfg_idx].split()
    # replace the first token that endswith .xyz
    replaced = False
    for k,tok in enumerate(cfg_line):
        if tok.lower().endswith(".xyz"):
            cfg_line[k] = start_xyz_name
            replaced = True
            break
    if not replaced:
        cfg_line.append(start_xyz_name)
    lines[cfg_idx] = " ".join(cfg_line)
    cass_inp_path.write_text("\n".join(lines) + "\n")
    print(f"[{now()}] Cassanda input file updated: box {old_val} -> {lines[bi+3]}, start xyz -> {start_xyz_name}")
    

    #optional backup_file
    if backup:
        backup_file(cass_inp_path)

def parse_molecule_files_from_cass_inp(cass_inp_path: Path):
    """Return list of tuples (filename, count) from '# Molecule_Files' block in Cassanda input"""
    text = cass_inp_path.read_text()
    lines = text.splitlines()
    start = None
    for i,ln in enumerate(lines):
        if ln.strip().startswith("# Molecule_Files"):
            start = i + 1
            break
    if start is None:
        return []
    mols = []
    for j in range(start, len(lines)):
        ln = lines[j].strip()
        # skip empty lines or comment lines starting with # or !
        if ln == "" or ln.startswith("#") or ln.startswith("!"):
            break
        parts = ln.split()
        fname = parts[0]
        cnt = None
        if len(parts) >= 2:
            try:
                cnt = int(float(parts[1]))
            except:
                cnt = None
        mols.append((fname, cnt))
    return mols

def parse_prp_for_last_counts(prp_path: Path, mol_files, skip_cols=2):
    """Parse .prp; return mapping molebasename -> last integer count (or None if not parsed)"""
    content = prp_path.read_text(errors='ignore')
    lines = content.splitlines()
    data_lines = []
    for ln in lines:
        if ln.strip().startswith("#"):
            continue
        if ln.strip() == "":
            continue
        data_lines.append(ln.strip())
    if not data_lines:
        raise RuntimeError(f"{prp_path} contains no data lines")
    last = data_lines[-1]
    parts = last.split()
    results = {}
    names = [Path(x[0]).stem for x in mol_files]
    for i, name in enumerate(names):
        col = skip_cols + i
        if col < len(parts):
            raw = parts[col]
            try:
                val = int(float(raw))
            except:
                try:
                    val = round(float(raw))
                    val = int(val)
                except:
                    val = None
            results[name] = val
        else:
            results[name] = None
    return results

def edit_xyz_to_gro_script(script_path: Path, update_counts: dict, box_nm: float, output_gro_name: str, backup=False):
    """Edit xyz_to_gro.py in-place to update molecules dict counts, box_size, and output_gro name.
       This attempts to find a 'molecules = {' block and replace values for matching keys.
    """
    text = script_path.read_text()
    orig = text

    # Replace molecules block - naive regex approach
    mblock_re = re.compile(r"(molecules\s*=\s*\{)(.*?)(\})", re.S)
    m = mblock_re.search(text)
    if m:
        inner = m.group(2)
        lines = inner.splitlines()
        new_lines = []
        for ln in lines:
            s = ln.strip()
            if s == "":
                new_lines.append(ln)
                continue
            keymatch = re.match(r"""['"]?([A-Za-z0-9_+-]+)['"]?\s*:\s*\(\s*([0-9]+)\s*,\s*([0-9]+)\s*\)\s*,?""", s)
            if keymatch:
                key = keymatch.group(1)
                a = int(keymatch.group(2))
                b = int(keymatch.group(3))
                if key in update_counts:
                    nb = int(update_counts[key])
                    newln = re.sub(r"\([^\)]*\)", f"({a}, {nb})", ln)
                    new_lines.append(newln)
                else:
                    new_lines.append(ln)
            else:
                new_lines.append(ln)
        new_inner = "\n".join(new_lines)
        text = text[:m.start(2)] + new_inner + text[m.end(2):]

    # Replace box_size
    box_re = re.compile(r"(\bbox_size\s*=\s*)([0-9]*\.?[0-9]+)")
    if box_re.search(text):
        text = box_re.sub(lambda mo: f"{mo.group(1)}{box_nm:.5f}", text, count=1)

    # Replace output_gro var
    out_re = re.compile(r"(\boutput_gro\s*=\s*)['\"].*?['\"]")
    if out_re.search(text):
        text = out_re.sub(lambda mo: f'{mo.group(1)}"{output_gro_name}"', text, count=1)

    if text != orig:
        #backup_file(script_path)
        script_path.write_text(text)
        print(f"[{now()}] Edited {script_path.name}: updated molecule counts, box size, and output_gro.")
        return True
    else:
        print(f"[{now()}] No changes made to {script_path.name}.")
        return False
        
    #optional backup_file
    if backup:
        backup_file(script_path)

def update_system_top(top_path: Path, new_counts: dict, backup=False):
    """Edit [ molecules ] block in system.top to set counts for new_counts (keys matched case-insensitive)"""
    text = top_path.read_text()
    lines = text.splitlines()
    start = None
    for i,ln in enumerate(lines):
        if ln.strip().lower().startswith("[ molecules ]"):
            start = i
            break
    if start is None:
        raise RuntimeError("No [ molecules ] section in system.top")
    out = lines[:start+1]
    j = start + 1
    while j < len(lines):
        ln = lines[j]
        if ln.strip().startswith("["):
            break
        if ln.strip() == "" or ln.strip().startswith(";"):
            out.append(ln)
            j += 1
            continue
        parts = ln.split()
        name = parts[0]
        matched_key = None
        for k in new_counts.keys():
            if k.lower() == name.lower():
                matched_key = k
                break
        if matched_key:
            out.append(f"{name}    {int(new_counts[matched_key])}")
        else:
            out.append(ln)
        j += 1
    # append the rest
    out.extend(lines[j:])
    #backup_file(top_path)
    top_path.write_text("\n".join(out) + "\n")
    print(f"[{now()}] system.top updated for molecule counts.")
    #def update_system_top(top_path, counts, backup=False):
    if backup:
        backup_file(top_path)

# ------------------------------
# High-level process functions
# ------------------------------
def apply_pbc_and_convert_gro_to_xyz(gmx_cmd, obabel_cmd, start_gro: Path, start_tpr: Path, dest_dir: Path, out_prefix: str):
    """Run trjconv PBC and obabel convert to xyz"""
    dest_pbc = dest_dir / f"{out_prefix}-pbc.gro"
    # gmx trjconv -f start_conf.gro -s start_conf.tpr -o start_conf-pbc.gro -pbc mol
    trjconv_cmd = f"echo 0 | {gmx_cmd} trjconv -f {start_gro} -s {start_tpr} -o {dest_pbc} -pbc mol"
    t0 = timenow()
    run_cmd(trjconv_cmd, cwd=dest_dir, shell=True)
    t1 = timenow()
    print(f"[{now()}] trjconv finished (elapsed {t1 - t0:.2f} s): {dest_pbc}")
    # obabel convert
    dest_xyz = dest_dir / f"{out_prefix}-pbc.xyz"
    run_cmd([obabel_cmd, "-i", "gro", str(dest_pbc), "-o", "xyz", "-O", str(dest_xyz)], cwd=dest_dir)
    t2 = timenow()
    print(f"[{now()}] obabel conversion finished (elapsed {t2 - t1:.2f} s): {dest_xyz}")
    return dest_pbc, dest_xyz, (t1 - t0, t2 - t1)

def run_cassandra_and_wait(cassandra_cmd, cass_inp_path: Path, cwd: Path, completion_string, check_interval=10, timeout=None, extra_failure_triggers=None):
    """Run Cassandra and wait until completion string appears in stdout or log.  Includes automatic runtime failure detection."""
    import re
    import time

    t0 = timenow()
    completion_norm = re.sub(r'\W+', '', completion_string.lower())
    found = False
    # ---------------------------------------------------------
    # BUILT-IN FAILURE TRIGGERS (engine-level protection)
    # ---------------------------------------------------------  
    #These strings are used to help the automation process detect where there is a failure, if any of these words are found, the cose will automatically fail.
    
    default_failure_triggers = [
        "fatal error",
        "error:",
        "stopping program",
        "cannot be computed",
        "segmentation fault",
        "forrtl",
        "nan",
    ]

    # Merge user-provided failure triggers (from YAML)
    user_triggers = [f.lower() for f in (extra_failure_triggers or [])]
    failure_triggers = default_failure_triggers + user_triggers          
    last_output = ""

    print(f"[{now()}] Launching Cassandra...")

    # start Cassandra process
    proc = subprocess.Popen([cassandra_cmd, str(cass_inp_path)],
                            cwd=cwd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            text=True)

    start_time = time.time()
    log_files = list(cwd.glob("*.log"))

    while True:
        # Here we monitor the Standard output - read any new stdout lines 
        if proc.stdout:
            for line in proc.stdout:
                last_output = line.strip()
                print(f"[CASS] {last_output}")
                norm_line = re.sub(r'\W+', '', last_output.lower())
               
                # Detect any of the success flags          
                if completion_norm in  norm_line:
                    found = True
                    break
                    
                    
                 # Detect any of the failure from the stdout       
                if any(err in last_output.lower() for err in failure_triggers):
                    proc.kill()
                    raise RuntimeError(f"Cassandra runtime error detected:\n{last_output}")                               
                
        if found:
            break

        #Here we check the log files for the success or failures too 
        
        for lg in log_files:
            try:
                txt = lg.read_text(errors='ignore').lower()
                
                # success
                
                if completion_norm in re.sub(r'\W+', '', txt):
                    found = True
                    break
                    
                for err in failure_triggers:
                    if err in txt:
                        proc.kill()
                        raise RuntimeError(
                            f"Cassandra error found in log file {lg.name}: '{err}'"
                        )                                        
                    
            except:
                continue

        if found:
            break
            
          # silent crash detection
            
        if proc.poll() is not None and not found:
            raise RuntimeError(
                f"Cassandra terminated unexpectedly.\nLast output:\n{last_output}"
            )      

        # check timeout
        if timeout and (time.time() - start_time) > timeout:
            proc.kill()
            raise RuntimeError(f"Cassandra did not finish within {timeout} seconds.")

        time.sleep(check_interval)

    proc.wait()  # ensure process fully exited
    t1 = timenow()
    print(f"[{now()}] Cassandra finished successfully (elapsed {t1 - t0:.2f}s).")
    return last_output, t1 - t0

def run_xyz_to_gro(python_exec, xyz_to_gro_script: Path, cwd: Path):
    """Run the user's xyz_to_gro.py inside cwd (script is copied to cwd prior)"""
    t0 = timenow()
    run_cmd([python_exec, str(xyz_to_gro_script)], cwd=cwd)
    t1 = timenow()
    print(f"[{now()}] xyz_to_gro.py finished (elapsed {t1 - t0:.2f} s).")
    return t1 - t0

def run_gromacs_min_md(gmx_cmd, min_mdp: Path, md_mdp: Path, input_gro: Path, top_file: Path, out_prefix: str, cwd: Path, extra_args=None):
    """Run grompp+mdrun for minimization and NPT. Returns path to final gro and tpr"""
    t0 = timenow()
    # gromacs energy minimization via  grompp em
    run_cmd([gmx_cmd, "grompp", "-f", str(min_mdp), "-c", str(input_gro), "-p", str(top_file), "-o", str(cwd / "em"), "-maxwarn", "10"], cwd=cwd)
    run_cmd([gmx_cmd, "mdrun", "-deffnm", "em"], cwd=cwd)

    # gromacs md run via grompp md
    run_cmd([gmx_cmd, "grompp", "-f", str(md_mdp), "-c", str(cwd / "em.gro"), "-p", str(top_file), "-o", str(cwd / out_prefix), "-maxwarn", "10"], cwd=cwd)

    # default GPU/CPU flags
    mdrun_cmd = [
        gmx_cmd, "mdrun", "-deffnm", out_prefix,
        "-pin", "on",
        "-nb", "gpu",
        "-ntmpi", "1",
        "-bonded", "gpu",
        "-pme", "gpu"
    ]

    # allow YAML to override/add -ntomp
    if extra_args:
        # if user provides ["-ntomp", "16"], it will be inserted correctly
        mdrun_cmd += extra_args

    run_cmd(mdrun_cmd, cwd=cwd)
    t1 = timenow()

    final_gro = cwd / f"{out_prefix}.gro"
    if not final_gro.exists():
        # attempt confout.gro fallback
        if (cwd / "confout.gro").exists():
            shutil.move(str(cwd / "confout.gro"), str(final_gro))
        else:
            raise RuntimeError("GROMACS did not produce final gro after relaxation")
    print(f"[{now()}] GROMACS min+md finished (elapsed {t1 - t0:.2f} s). final {final_gro}")
    return final_gro, (t1 - t0)

def trjconv_pbc_and_obabel_convert(gmx_cmd, obabel_cmd, tpr_path: Path, gro_path: Path, out_prefix: str, cwd: Path):
    """Run echo 0 | gmx trjconv -f gro -s tpr -o out-pbc.gro -pbc mol and then obabel"""
    t0 = timenow()
    cmd = f"echo 0 | {gmx_cmd} trjconv -f {gro_path} -s {tpr_path} -o {cwd / (out_prefix + '-pbc.gro')} -pbc mol"
    run_cmd(cmd, cwd=cwd, shell=True)
    # obabel
    run_cmd([obabel_cmd, "-i", "gro", str(cwd / (out_prefix + "-pbc.gro")), "-o", "xyz", "-O", str(cwd / (out_prefix + "-pbc.xyz"))], cwd=cwd)
    t1 = timenow()
    print(f"[{now()}] trjconv pbc + obabel done (elapsed {t1 - t0:.2f} s).")
    return cwd / (out_prefix + "-pbc.gro"), cwd / (out_prefix + "-pbc.xyz"), (t1 - t0)
    

def detect_prp_skip_columns(cass_inp_path: Path):
    """
    Detect the Nmols column index from Cassandra Property_Info block.

    Returns:
        skip_cols (int):
            Number of columns before Nmols in .prp file.
    """

    text = cass_inp_path.read_text()
    lines = text.splitlines()

    start = None
    
    # Locate Property_Info block
    for i, ln in enumerate(lines):
        if ln.strip().lower().startswith("# property_info"):
            start = i + 1
            break

    if start is None:
        raise RuntimeError(
            "Could not locate '# Property_Info' block in Cassandra input."
        )

    props = []

    for j in range(start, len(lines)):
        ln = lines[j].strip()
        
        # next block reached
        if ln.startswith("#"):
            break
            
        # skip empty lines safely
        if ln:
            props.append(ln.lower())

    if not props:
        raise RuntimeError(
            "Property_Info block found but contains no properties."
        )

    try:
        nmols_index = next(
            i for i, p in enumerate(props)
            if p in ("nmol", "nmols")
        )
    except StopIteration:
        raise RuntimeError(
            "Nmol/Nmols property not found in Property_Info block."
        )
        
     # STEP column always exists in .prp
    skip_cols = 1 + nmols_index

    print(f"[{now()}] Cassandra Property_Info detected:")
    print(f"[{now()}] Properties: {props}")
    print(f"[{now()}] Nmols position: {nmols_index}")
    print(f"[{now()}] Auto-detected skip columns = {skip_cols}")

    return skip_cols

# ------------------------------
# Utilities for restart
# ------------------------------
def cycle_completed(cyc_dir: Path, cyc: int):
    """Check if a cycle was fully completed by verifying final gro/tpr files exist."""
    tag = f"cycle{cyc:03d}"
    gro_file = cyc_dir / f"{tag}.gro"
    tpr_file = cyc_dir / f"{tag}.tpr"
    return gro_file.exists() and tpr_file.exists()


def detect_last_completed_cycle(work_dir: Path):
    """Scan cycle directories to find the last fully completed cycle.
       If a cycle exists but is incomplete, remove it for safety."""
    last_ok = 0
    for d in sorted(work_dir.glob("cycle_*")):
        cyc = int(d.name.split("_")[1])
        if cycle_completed(d, cyc):
            last_ok = cyc
        else:
            print(f"[WARN] Cycle {cyc} incomplete - removing directory for safety")
            shutil.rmtree(d)
            break
    return last_ok
    

# ------------------------------
#Banner function
# ------------------------------

def print_banner():

    banner = r"""
 _    _       _____ _           __   __
| |  | |     / ____(_)          \ \ / /
| |__| |_   _| (___  _ _ __ ___  \ V / 
|  __  | | | |\___ \| | '_ ` _ \  > <  
| |  | | |_| |____) | | | | | | |/ . \ 
|_|  |_|\__, |_____/|_|_| |_| |_/_/ \_\
         __/ |                         
        |___/                          

   HySimX : A Hybrid MD-GCMC Simulation Framework
            for eXtended Adsorption and Structural Relaxation in Host-Guest Systems
"""

    print("\n" + "="*80)
    print(banner)
    print("="*80 + "\n")
    print(" Author  : Ifeanyi Samuel Ibeto")
    print(" Version : 1.0")
    print(" Date    : November 2025")
    print("-"*80)
    print(" Engines : GROMACS (MD)  |  Cassandra (GCMC)")
    print(" Tools   : Python  |  OpenBabel")
    print("-"*80)
    print(" MD says : 'Lets relax the host structure.'")
    print(" MC says : 'Sample the guest molecules population.'")
    print(" HySimX  : 'Cycle them until equilibrium.'")    
    print("="*80 + "\n")       
    print("If you find this framework useful please Cite: Ibeto et al., HySimX Framework (2026)")
    print("-"*80)
    print("Initializing HySimX...")
    print("-"*80)
    print()

  

# ------------------------------
# Main cycle manager
# ------------------------------


# Print banner at the very start
print_banner()



def main():
    import argparse    
    import logging
# ------------------------------
    # Parse command-line args
    # ------------------------------
    parser = argparse.ArgumentParser(description="MD-MC Hybrid Cycle Manager")
    parser.add_argument("config.yml", help="YAML configurational file path containing all the run details")
    parser.add_argument("--restart", nargs="?", type=int, const=None,
                        help="Optional: number of additional cycles to run beyond YAML cycles. "
                             "If no number provided, continue until YAML cycles are met.")
    args = parser.parse_args()
      
        
    cfg_path = Path(sys.argv[1])
    cfg = load_config(cfg_path)
    backup_enabled = cfg.get('backup', False)  #optional backup-not really necessary
    work_dir = Path.cwd().resolve()
    cycles = int(cfg.get('cycles', 1))
    yaml_cycles = int(cfg.get('cycles', 1))
    
    # ------------------------------
    # Detect last fully completed cycle from previous runs
    # ------------------------------
    last_done = detect_last_completed_cycle(work_dir)
    
    
    # ------------------------------
    # Load existing checkpoint (append-safe)
    # ------------------------------
    checkpoint_file = work_dir / "checkpoint.json"
    if checkpoint_file.exists():
        try:
            checkpoint = json.loads(checkpoint_file.read_text())
            if "history" not in checkpoint:
                checkpoint["history"] = []
            print(f"[{now()}] Existing checkpoint loaded. Last cycle: {last_done}")
        except Exception as e:
            print(f"[{now()}] WARNING: Could not read checkpoint: {e}. Starting fresh.")
            checkpoint = {"last_cycle": 0, "history": []}
    else:
        checkpoint = {"last_cycle": 0, "history": []}
    
    
    

    # ------------------------------
    # Determine target cycles based on restart flag
    # ------------------------------
    if args.restart is not None:
        # --restart provided
        if args.restart is None:
            # User did not specify a number: run until YAML cycles
            target_cycles = max(yaml_cycles, last_done)
            extra_cycles = target_cycles - last_done
        else:
            extra_cycles = int(args.restart)
            target_cycles = last_done + extra_cycles
    else:
        # Normal run
        extra_cycles = 0
        target_cycles = yaml_cycles

    # ------------------------------
    # STARTUP MESSAGES: Inform user
    # ------------------------------
#    print("\n" + "=" * 70)
#    print(f"[{now()}] MD-MC Hybrid Cycle Manager Starting")
#    print("=" * 70)
    print(f"[{now()}] Loaded config from: {cfg_path}")
    print(f"[{now()}] Work directory: {work_dir}")
    print(f"[{now()}] YAML-specified cycles: {yaml_cycles}")
    print(f"[{now()}] Last fully completed cycle detected: {last_done}")

    if last_done == 0:
        print(f"[{now()}] Starting a fresh run from cycle 1 to {target_cycles}")
    elif last_done >= yaml_cycles and args.restart is None:
        print(f"[{now()}] All YAML-specified cycles already completed. No action taken.")
        return
    elif args.restart is None:
        remaining = yaml_cycles - last_done
        print(f"[{now()}] Continuing run: {remaining} remaining cycles to reach YAML target ({yaml_cycles})")
    else:
        print(f"[{now()}] Restart detected. Resuming from cycle {last_done + 1} to {target_cycles}")

    print(f"[{now()}] Total cycles to execute in this run: {target_cycles - last_done}")
    print("=" * 70 + "\n")


    # ------------------------------
    # Determine start cycle
    # ------------------------------
 
    start_cycle = last_done + 1
 #   print(f"[{now()}] Resuming from cycle {start_cycle} to {target_cycles}")

    # ------------------------------
    # validate starting files
    # ------------------------------

    print(f"[{now()}] Validating starting files...")
    mol_files = validate_start_files(cfg, work_dir)
    print(f"[{now()}] Starting files validated.")

    mapping = cfg['molecules']['mapping']
    gases = cfg['molecules']['gases']
    adsorbents = cfg['molecules']['adsorbent']

    checkpoint = {"last_cycle": 0, "history": []}
    
    # ------------------------------
    # Detect last completed cycle
    # ------------------------------
    last_done = detect_last_completed_cycle(work_dir)
    if last_done >= target_cycles:
        print(f"[{now()}] All {target_cycles} cycles already complete. Nothing to do.")
        return

    start_cycle = last_done + 1
    print(f"[{now()}] Resuming from cycle {start_cycle}")
           
    

    # ------------------------------
    # LOOP OVER CYCLES
    # ------------------------------
    for cyc in range(start_cycle, target_cycles + 1):
        print("\n" + "="*60)
        print("\n" + "="*60)
        #print(f"[{now()}] <<<<<<<<<< ENTERING CYCLE {cyc}  >>>>>>>>>>")
        print(f"           <<<<<<<<<<  CYCLE {cyc}  >>>>>>>>>>")
        print("\n" + "="*60)
        print("\n" + "="*60)
        
        cycle_start = timenow()

        cyc_dir = work_dir / f"cycle_{cyc:03d}"
        cyc_dir.mkdir(parents=True, exist_ok=True)
#        print(f"[{now()}] Created cycle directory: {cyc_dir}")

        # copy ALL.inp and mcf + species
#        print(f"[{now()}] Copying ALL.inp, mcf files, and species folders...")
        cass_dir = (work_dir / cfg['cassandra']['cass_inp']).parent

        cass_inp_name = Path(cfg['cassandra']['cass_inp']).name
        shutil.copy2(work_dir / cfg['cassandra']['cass_inp'], cyc_dir / cass_inp_name)
        
        
        for fname, cnt in mol_files:
            src1 = cass_dir / fname
            src2 = work_dir / fname
            if src1.exists():
#                print(f"[{now()}] Copying {src1} -> {cyc_dir}")
                shutil.copy2(src1, cyc_dir / fname)
            elif src2.exists():
#                print(f"[{now()}] Copying {src2} -> {cyc_dir}")
                shutil.copy2(src2, cyc_dir / fname)

        species_prefix = cfg['cassandra'].get('species_folder_prefix', 'species')
        for sp in cass_dir.glob(f"{species_prefix}*"):
            if sp.is_dir():
#                print(f"[{now()}] Copying species folder {sp.name}")
                dest = cyc_dir / sp.name
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(sp, dest)

        # copy scripts + topology + mdp
#        print(f"[{now()}] Copying xyz_to_gro.py, topology, and mdp files...")
        for name in [cfg['xyz_to_gro']['script'], cfg['gromacs']['topology'], cfg['gromacs']['min_mdp'], cfg['gromacs']['md_mdp']]:
            src = work_dir / name
            if src.exists():
#                print(f"[{now()}] Copying {src} -> {cyc_dir}")
                shutil.copy2(src, cyc_dir / src.name)

        # --------------------------------------------------------------
        # STEP 1: Choose starting gro/tpr for THIS cycle
        # --------------------------------------------------------------
#        print(f"[{now()}] Selecting start files for cycle {cyc}...")

        if cyc == 1:
#            print(f"[{now()}] Cycle 1 -> using initial start files.")
            start_gro = work_dir / cfg['start']['gro']
            start_tpr = work_dir / cfg['start']['tpr']
            out_prefix = "start_conf"
        else:
#            print(f"[{now()}] Cycle {cyc} -> using previous cycle outputs.")
            prev_tag = f"cycle{(cyc-1):03d}"
            prev_dir = work_dir / f"cycle_{cyc-1:03d}"
            prev_gro = prev_dir / f"{prev_tag}.gro"
            prev_tpr = prev_dir / f"{prev_tag}.tpr"

#            print(f"[{now()}] Checking: {prev_gro}, {prev_tpr}")
            if not prev_gro.exists() or not prev_tpr.exists():
                raise RuntimeError(
                    f"Missing previous cycle files for cycle {cyc}: {prev_gro} or {prev_tpr}"
                )

            start_gro = prev_gro
            start_tpr = prev_tpr
            out_prefix = f"cycle{cyc-1:03d}"   #here the new file is renamed for pbc application and conversion to xyz

#        print(f"[{now()}] Running PBC + xyz conversion...")
        pbc_gro, pbc_xyz, times = apply_pbc_and_convert_gro_to_xyz(
            cfg['executables']['gmx'],
            cfg['executables']['obabel'],
            start_gro, start_tpr,
            cyc_dir, out_prefix
        )
#        print(f"[{now()}] PBC gro: {pbc_gro}")
#        print(f"[{now()}] PBC xyz: {pbc_xyz}")

        # --------------------------------------------------------------
        # STEP 2: update Cassandra input file 
        # --------------------------------------------------------------
#        print(f"[{now()}] Updating Cassandra input file ...")
        box_nm = gro_box_size_nm(pbc_gro)
        box_A = box_nm * 10.0
        
       # backup_file(cyc_dir / "ALL.inp")
        cass_inp_name = Path(cfg['cassandra']['cass_inp']).name
        cass_inp_path = cyc_dir / cass_inp_name
        
        if cfg.get('backup', False):  
            backup_file(cass_inp_path)
        
                      
        edit_cassandra_inp_replace_box_and_start(cass_inp_path, box_A, pbc_xyz.name,backup=backup_enabled)
#        print(f"[{now()}] ALL.inp updated with box size {box_A} A")

        # --------------------------------------------------------------
        # STEP 3: run Cassandra
        # --------------------------------------------------------------
        print(f"[{now()}] Running Cassandra...")
        cass_out, cass_time = run_cassandra_and_wait(
            cfg['executables']['cassandra'],
            cass_inp_path,
            cyc_dir,
            cfg['cassandra']['completion_string'],
            extra_failure_triggers=cfg['cassandra'].get('extra_failure_triggers', [])
        )
        print(f"[{now()}] Cassandra finished in {cass_time:.2f} s")

        # --------------------------------------------------------------
        # STEP 4: parse prp
        # --------------------------------------------------------------
#        print(f"[{now()}] Parsing .prp file...")
        prp_candidates = list(cyc_dir.glob(cfg['cassandra']['prp_pattern']))
        if not prp_candidates:
            raise RuntimeError("No .prp file found after Cassandra run in " + str(cyc_dir))

        prp = sorted(prp_candidates)[-1]
#        print(f"[{now()}] Found {prp.name}")

        
        # Automatically detect Nmols column position
        skip_cols = detect_prp_skip_columns(cass_inp_path)
        
        counts = parse_prp_for_last_counts(
            prp,
            mol_files,
            skip_cols=skip_cols
        )
                               
        #counts = parse_prp_for_last_counts(prp, mol_files, skip_cols=cfg.get('tolerances', {}).get('skip_prp_columns', 2))
#        print(f"[{now()}] Parsed counts: {counts}")
        
        gas_updates = {}
        for g in gases:
            if g in counts and counts[g] is not None:
                gas_updates[g] = counts[g]
            else:
                raise RuntimeError(f"Gas '{g}' count missing from prp")

        # --------------------------------------------------------------
        # STEP 5: update xyz_to_gro.py
        # --------------------------------------------------------------
#        print(f"[{now()}] Updating xyz_to_gro.py for new counts...")
        script_name = Path(cfg['xyz_to_gro']['script']).name
        script_in_cycle = cyc_dir / script_name
        output_gro_name = f"output-{cyc:03d}.gro"
        edit_xyz_to_gro_script(script_in_cycle, gas_updates, box_nm, output_gro_name,backup=backup_enabled)
#        print(f"[{now()}] xyz_to_gro.py updated.")

        # --------------------------------------------------------------
        # STEP 6: run xyz_to_gro.py
        # --------------------------------------------------------------
#        print(f"[{now()}] Running xyz_to_gro.py...")
        run_xyz_to_gro(cfg['executables']['python'], script_in_cycle, cwd=cyc_dir)

        out_gro = cyc_dir / output_gro_name
        if not out_gro.exists():
            print(f"[{now()}] WARNING: expected gro not found, searching...")
            all_gros = list(cyc_dir.glob("*.gro"))
            if not all_gros:
                raise RuntimeError("xyz_to_gro.py did not produce any gro file.")
            out_gro = sorted(all_gros)[-1]

#        print(f"[{now()}] xyz_to_gro output gro: {out_gro}")

        # STEP 7: update topology
#        print(f"[{now()}] Updating system.top...")
        top_in_cycle = cyc_dir / Path(cfg['gromacs']['topology']).name

        new_counts_map = {k: v for k, v in counts.items() if v is not None}
        update_system_top(top_in_cycle, new_counts_map,backup=backup_enabled)
#        print(f"[{now()}] system.top updated.")

        # --------------------------------------------------------------
        # STEP 8: run GROMACS min+npt
        # --------------------------------------------------------------
        print(f"[{now()}] Running GROMACS minimization + NPT...")
        extra = cfg['gromacs'].get('md_extra_args', [])
        final_gro, md_time = run_gromacs_min_md(
            cfg['executables']['gmx'],
            cyc_dir / Path(cfg['gromacs']['min_mdp']).name,
            cyc_dir / Path(cfg['gromacs']['md_mdp']).name,
            out_gro, top_in_cycle,
            f"cycle{cyc:03d}",
            cyc_dir,
            extra_args=extra
        )
        print(f"[{now()}] GROMACS finished in {md_time:.2f} s")
                       

        # --------------------------------------------------------------
        # RECORD CHECKPOINT (append-safe)
        # --------------------------------------------------------------
        cycle_end = timenow()
        cycle_elapsed = cycle_end - cycle_start

        # Append new cycle entry to history
        checkpoint_entry = {
            "cycle": cyc,
            "start_time": now(),
            "elapsed_s": cycle_elapsed,
            "counts": counts,
            "pbc_xyz": str(pbc_xyz),
            "pbc_gro": str(pbc_gro)
        }
        checkpoint["history"].append(checkpoint_entry)

        # Update last_cycle field to current cycle
        checkpoint["last_cycle"] = cyc

        # Write checkpoint back safely (preserve all history)
        checkpoint_file.write_text(json.dumps(checkpoint, indent=2))
        
        
        
        print("\n" + "="*60)
        print(f"[{now()}] <<< LEAVING CYCLE {cyc} - elapsed {cycle_elapsed:.2f} s")
        print("\n" + "="*60)
        
    print(f"[{now()}] All cycles finished. Last cycle = {checkpoint.get('last_cycle')}")

if __name__ == "__main__":
    main()
        

