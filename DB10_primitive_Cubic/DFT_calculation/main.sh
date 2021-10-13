#!/bin/bash
# The program is used to generate the DFT database of primitive cell.
# lei.zhang@rug.nl, 9th June 2021.

#--------------INITIAL SETTINGS------------
# Specify the file that contains atomic coordinates
atom_coord="lattice"
# Set k-spacing (unit: 1/(2*pi))
kspacing=0.015
# Specify kinetic energy and electron density cutoff
ECUTWFC='90'
ECUTRHO='1080'
# Prefix of the folder name
PREFIX='Fe_prim'

# Get the number of total data points
Tot=`wc -l <  ${atom_coord}`

# -----------------MAIN LOOP--------------------
for i in $(seq 0 1 $Tot)
do
# Define file names for QE
QE_INPUT="$PREFIX.$i.in"
QE_OUTPUT="$PREFIX.$i.out"
submitfile="submit_prim.$i"

# ---------------READ DATA-----------------
# Read cell parameters
a1=$(awk 'NR=='$i' {print $1}' ${atom_coord})
a2=$(awk 'NR=='$i' {print $2}' ${atom_coord})
a3=$(awk 'NR=='$i' {print $3}' ${atom_coord})
b1=$(awk 'NR=='$i' {print $4}' ${atom_coord})
b2=$(awk 'NR=='$i' {print $5}' ${atom_coord})
b3=$(awk 'NR=='$i' {print $6}' ${atom_coord})
c1=$(awk 'NR=='$i' {print $7}' ${atom_coord})
c2=$(awk 'NR=='$i' {print $8}' ${atom_coord})
c3=$(awk 'NR=='$i' {print $9}' ${atom_coord})

# Calculate the k points density
lengtha=`awk "BEGIN {x=$a1*$a1+$a2*$a2+$a3*$a3; y=sqrt(x); z=1/y; mu=z/$kspacing; print mu}"`
lengthb=`awk "BEGIN {x=$b1*$b1+$b2*$b2+$b3*$b3; y=sqrt(x); z=1/y; mu=z/$kspacing; print mu}"`
lengthc=`awk "BEGIN {x=$c1*$c1+$c2*$c2+$c3*$c3; y=sqrt(x); z=1/y; mu=z/$kspacing; print mu}"`
k1=`echo $lengtha | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}'`
k2=`echo $lengthb | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}'`
k3=`echo $lengthc | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}'`

#------------------------
# Create the input file for QE
#--------------------------------------------
cat > $QE_INPUT << EOF
&control
    calculation = 'scf',
    verbosity = 'high'
    restart_mode= 'from_scratch'
    prefix = 'Fe_prime'
    outdir = './out'
    pseudo_dir = './pseudo'
    tprnfor =.true.
    tstress =.true.
    disk_io='none'
    wf_collect = .false.
    max_seconds= 82800
 /
&system
    ibrav = 0
    nat = 1
    ntyp =1
    ecutwfc = $ECUTWFC,
    ecutrho = $ECUTRHO,
    occupations = 'smearing',
    smearing = 'marzari-vanderbilt',
    degauss = 0.01
    nspin =2
    starting_magnetization(1)=0.4
 /

 &electrons
 electron_maxstep = 300,
 conv_thr=1e-10,
 mixing_beta=0.2,
 /

ATOMIC_SPECIES
Fe 55.845 Fe.pbe-spn-rrkjus_psl.0.2.1.UPF

ATOMIC_POSITIONS angstrom
Fe 0.000000 0.000000 0.000000

K_POINTS automatic
$k1 $k2 $k3 1 1 1

CELL_PARAMETERS angstrom
$a1 $a2 $a3
$b1 $b2 $b3
$c1 $c2 $c3
EOF
#--------------------------------------------
# Create the file for submitting jobs by using batch
#--------------------------------------------
cat >${submitfile} <<EOF
#!/bin/bash
#SBATCH --job-name=PM_$i
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=30:00:00

module load gcc
module load intel/compiler/64/2020/19.1.1
module load intel/mpi/64/2018/4.274
module load intel/mkl/64/2018/4.274
module load python/python-3.6.0
module load qe/6.2

mpirun -np 32 pw.x -input ${QE_INPUT} > ${QE_OUTPUT}

EOF
#--------------------------------------------

# Move all file to a new dir
mkdir prim.$i
mv ${QE_INPUT} prim.$i
mv ${submitfile} prim.$i
# copy a pseudopotential file to target folder
cp -r pseudo prim.$i
# Sumbit the job
cd prim.$i
sbatch ${submitfile}
cd ..

echo 'Job submitted successfully.'
echo 'Job type: primitive cell calculation'
echo 'job id:' $i
sleep 1

done
