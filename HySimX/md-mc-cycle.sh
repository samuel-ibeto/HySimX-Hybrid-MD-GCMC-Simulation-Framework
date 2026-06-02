#!/bin/bash
#SBATCH --qos main
#SBATCH -J 2C-co-cy
#SBATCH --mail-user=isibeto@crimson.ua.edu
#SBATCH --mail-type=END
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --nodelist=gpu06

# ---- Load required modules ----
if ! shopt -q login_shell; then
  if [ -d /etc/profile.d ]; then
          for i in /etc/profile.d/*.sh; do
                  if [ -r $i ]; then
                          . $i
                  fi
          done
  fi
fi

module use /data3/share/modules

module load gromacs/2022.4



# ---- Activate your Python environment ----
#source ~/.bashrc
#conda activate mdmc     # OR replace with your environment name

# ---- Navigate to project directory ----
#cd /data3/samuel/PROJECTS/GAS-ADSORPTION/MD-MC-MD-framework/cycle-test/      # <-- EDIT THIS PATH

# ---- Run the MD–MC cycle driver ----

#python md-mc-cycle.py config.yml --restart 1

python -u  main.py config.yml
