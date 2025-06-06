#Gromacs preparation and result management

BASE_PDB = test
SOLV_PDB = solv_${BASE_PDB}
IONS_PDB = ions_${BASE_PDB}
INFILE = eq
RUNFILE=step5_1

#convert PDB to gro and prepare files
setup:
		#create gromacs representation
		gmx pdb2gmx -f ${BASE_PDB}.pdb -o ${BASE_PDB}.gro -i ${BASE_PDB}.itp -ff charmm36m -water tip3p -p ${BASE_PDB}.top
		#set system boundaries
		gmx editconf -f ${BASE_PDB}.gro -c -d 1.0 -bt triclinic -o out
		# Do not forget to include any needed forcefields into the topology file: ${BASE_PDB}.top
		kate ${BASE_PDB}.top
ions:
		#make iond.mdp
		echo "define                  = -DPOSRES " > ions.mdp
		echo "integrator              = steep" >> ions.mdp
		echo "emtol                   = 500.0" >> ions.mdp
		echo "nsteps                  = 50000" >> ions.mdp
		echo "nstlist                 = 10" >> ions.mdp
		echo "cutoff-scheme           = Verlet" >> ions.mdp
		echo "rlist                   = 1.2" >> ions.mdp
		echo "vdwtype                 = Cut-off" >> ions.mdp
		echo "vdw-modifier            = Force-switch" >> ions.mdp
		echo "rvdw_switch             = 1.0" >> ions.mdp
		echo "rvdw                    = 1.2" >> ions.mdp
		echo "coulombtype             = PME" >> ions.mdp
		echo "rcoulomb                = 1.2" >> ions.mdp
		echo "constraints             = h-bonds" >> ions.mdp
		echo "constraint_algorithm    = LINCS" >> ions.mdp
		# prepare solvated system, needs minimization
		#solvate
		cp ${BASE_PDB}.top topol.top
		gmx solvate -cp out.gro -cs spc216.gro -o ${SOLV_PDB}.gro -p topol.top
		#add ions
		gmx grompp -f ions.mdp -c ${SOLV_PDB}.gro -p topol.top -o ions.tpr -r ${SOLV_PDB}.gro -maxwarn 2
		#add 0.15 M KCl
		gmx genion -s ions.tpr -o ${IONS_PDB}.gro -p topol.top -pname K -nname CL -neutral -conc 0.15
index:
		#Make an index file. Select protein | ligand. Press enter and q
		gmx make_ndx -f ${IONS_PDB}.gro -o index.ndx

mini:
		#make emin.mdp
		echo "define                  = -DPOSRES " > emin.mdp
		echo "integrator              = steep" >> emin.mdp
		echo "emtol                   = 500.0" >> emin.mdp
		echo "nsteps                  = 50000" >> emin.mdp
		echo "nstlist                 = 10" >> emin.mdp
		echo "cutoff-scheme           = Verlet" >> emin.mdp
		echo "rlist                   = 1.2" >> emin.mdp
		echo "vdwtype                 = Cut-off" >> emin.mdp
		echo "vdw-modifier            = Force-switch" >> emin.mdp
		echo "rvdw_switch             = 1.0" >> emin.mdp
		echo "rvdw                    = 1.2" >> emin.mdp
		echo "coulombtype             = PME" >> emin.mdp
		echo "rcoulomb                = 1.2" >> emin.mdp
		echo "constraints             = h-bonds" >> emin.mdp
		echo "constraint_algorithm    = LINCS" >> emin.mdp
		gmx grompp -f emin.mdp -o emin.tpr -c ${IONS_PDB}.gro -r ${IONS_PDB}.gro -p topol.top -maxwarn 2 -n index.ndx
		gmx mdrun -nt 10 -pin on -v -deffnm emin &&    echo 'Potential' | gmx energy -f emin.edr -o emin_potential.xvg

showmin:
		xmgrace emin_potential.xvg

eq:
		#1_equilibration of the system
		#sed 's/^/echo "/' e.txt | sed 's/$/ ">> eq.mdp/'
		echo "define                  = -DPOSRES "> eq.mdp
		echo "integrator              = md ">> eq.mdp
		echo "dt                      = 0.001 ">> eq.mdp
		echo "nsteps                  = 125000 ">> eq.mdp
		echo "nstxout-compressed      = 5000 ">> eq.mdp
		echo "nstxout                 = 0 ">> eq.mdp
		echo "nstvout                 = 0 ">> eq.mdp
		echo "nstfout                 = 0 ">> eq.mdp
		echo "nstcalcenergy           = 100 ">> eq.mdp
		echo "nstenergy               = 1000 ">> eq.mdp
		echo "nstlog                  = 1000 ">> eq.mdp
		echo "; ">> eq.mdp
		echo "cutoff-scheme           = Verlet ">> eq.mdp
		echo "nstlist                 = 20 ">> eq.mdp
		echo "rlist                   = 1.2 ">> eq.mdp
		echo "vdwtype                 = Cut-off ">> eq.mdp
		echo "vdw-modifier            = Force-switch ">> eq.mdp
		echo "rvdw_switch             = 1.0 ">> eq.mdp
		echo "rvdw                    = 1.2 ">> eq.mdp
		echo "coulombtype             = PME ">> eq.mdp
		echo "rcoulomb                = 1.2 ">> eq.mdp
		echo "; ">> eq.mdp
		echo "tcoupl                  = v-rescale ">> eq.mdp
		echo "tc_grps                 = SYSTEM ">> eq.mdp
		echo "tau_t                   = 1.0 ">> eq.mdp
		echo "ref_t                   = 310 ">> eq.mdp
		echo "; ">> eq.mdp
		echo "constraints             = h-bonds ">> eq.mdp
		echo "constraint_algorithm    = LINCS ">> eq.mdp
		echo "; ">> eq.mdp
		echo "nstcomm                 = 100 ">> eq.mdp
		echo "comm_mode               = linear ">> eq.mdp
		echo "comm_grps               = SYSTEM ">> eq.mdp
		echo "; ">> eq.mdp
		echo "gen-vel                 = yes ">> eq.mdp
		echo "gen-temp                = 310 ">> eq.mdp
		echo "gen-seed                = -1 ">> eq.mdp
		gmx grompp -f eq.mdp -o eq.tpr -c emin.gro -r emin.gro -p topol.top -n index.ndx
		gmx mdrun -v -pin on -nt 10 -nb gpu -nstlist 400 -deffnm eq

checkeq:
		#energy temp stabilization
		gmx energy -f eq.edr -o temperature.xvg
		# pressure stabilization
		gmx energy -f eq.edr -o pressure.xvg

fixtc:
		#usage: make analysis RUNFILE=path_no_extension
		#fix PBC errors of trajectory xtc files
		#gmx trjconv -f ${RUNFILE}.xtc -s ${RUNFILE}.tpr -pbc mol -boxcenter zero -o temp1.xtc -n index.ndx
		#gmx trjconv -f ${RUNFILE}.xtc -s ${RUNFILE}.tpr -pbc whole -ur compact -center -o temp2.xtc -n index.ndx
		gmx trjconv -f ${RUNFILE}.xtc -s ${RUNFILE}.tpr -pbc nojump -boxcenter zero -o c_${RUNFILE}.xtc -n index.ndx



startpdb:
		#save first frame as a reference for later analysis
		gmx trjconv -s ${RUNFILE}.tpr -f ${RUNFILE}.xtc -o start.pdb -dump 0

analysis:
		#usage: make analysis INFILE=path_no_extension
		#RMSD
		gmx rms -s ${INFILE}.tpr -f ${INFILE}.xtc -o rmsd.xvg -tu ns
		#GYR
		gmx gyrate -s ${INFILE}.tpr -f ${INFILE}.xtc -o gyrate.xvg

combine:
		#concatenate all XTC files setting the correct time interactively
		gmx trjcat -f *.xtc -o comb.xtc -settime
