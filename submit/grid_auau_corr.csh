#!/bin/csh

#  A Slightly ugly way to make it as simple as possible to run different analyses,
#  And to make it easier to change variables and keep the book keeping simple
#  Command Line Arugments:
#  [1]: the analysis: either jet or dijet
#  [2]: if set to 'default', uses original A_j parameterization ( or jet-hadron )
#  [2]: if not using default settings, sets whether to use particle-by-particle
#       efficiency corrections
#
#  -----Di_jet selection criteria-----
#  [3]: whether to require a trigger in your leading jet
#  [4]: subleading jet minimum pt ( for jet hadron, set to zero )
#  [5]: leading jet minimum pt ( jet min pt for jet-hadron )
#  [6]: jet pt max ( global maximum, both dijet and jet )
#  [7]: jet radius for clustering algorithm
#  [8]: hard constituent pt cut
#  [8]: bins in Eta for correlation histograms
#  [9]: bins in phi for correlation histograms
#
#  Output names and locations are generated by the script, and correspond to the above variables

# first make sure program is updated and exists
 make bin/auau_correlation || exit

if ( $2 == '-h') then
echo 'parameters:'
echo 'for defaults use "(di)/jet default"'
echo '1: analysis type [dijet/jet] (default: dijet)'
echo '2: use tracking efficiency corrections [true/false] (default: true)'
echo '3: require trigger in leading jet [true/false] (default: true)'
echo '4: subleading jet min pt (default: 10)'
echo '5: leading jet min pt (default: 20)'
echo '6: jet pt max (default: 100)'
echo '7: jet resolution parameter (default: 0.4'
echo '8: bins in correlation histograms in eta (default: 22)'
echo '9: bins in correlation histograms in phi (default: 22)'
endif

set ExecPath = `pwd`
set analysis = $1
set execute = './bin/auau_correlation'
set base = /nfs/rhi/STAR/Data/CleanAuAuY7/Clean

if ( $# != "10" && !( $2 == 'default' ) ) then
	echo 'Error: illegal number of parameters (-h for help)'
	exit
endif

echo $analysis
if ( $analysis != 'dijet' && $analysis != 'jet' ) then
	echo 'Error: unknown analysis type'
	exit
endif

# Arguments
set useEfficiency = $2
set triggerCoincidence = $3
set subLeadPtMin = $4
set leadPtMin = $5
set jetPtMax = $6
set jetRadius = $7
set constPtCut = $8
set binsEta = $9
set binsPhi = $10

if ( $2 == 'default' ) then
	set useEfficiency = 'true'
	set triggerCoincidence = 'true'
	if ( $analysis == 'dijet' ) then
		set subLeadPtMin = 10.0
		set leadPtMin = 20.0
		set jetPtMax = 100.0
	else if ( $analysis == 'jet' ) then
		set subLeadPtMin = 0.0
		set leadPtMin = 20.0
		set jetPtMax = 100.0
	endif
  endif
	set jetRadius = 0.4
  set constPtCut = 2.0
  set binsEta = 22
  set binsPhi = 22
endif

# Create the folder name for output
set outFile = ${analysis}
set outFile = ${outFile}_trigger_${triggerCoincidence}_eff_${useEfficiency}_lead_${leadPtMin}_sub_${subLeadPtMin}_max_${jetPtMax}_rad_${jetRadius}_hardpt_${constPtCut}_eta_${binsEta}_phi_${binsPhi}
# Make the directories since they may not exist...
if ( ! -d out/${analysis}/${outFile} ) then
mkdir -p out/${analysis}/${outFile}
mkdir -p out/${analysis}/${outFile}/correlations
mkdir -p out/${analysis}/${outFile}/tree
mkdir -p out/${analysis}/${outFile}/mixing
endif

if ( ! -d log/${analysis}/${outFile} ) then
mkdir -p log/${analysis}/${outFile}
endif

# Now Submit jobs for each data file
foreach input ( ${base}* )

# Create the output file base name
set OutBase = `basename $input | sed 's/.root//g'`

# Make the output names and path
set outLocation = "out/${analysis}/${outFile}/"
set outName = correlations/corr_${OutBase}.root
set outNameTree = tree/tree_${OutBase}.root

# Input files
set Files = ${input}

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out
set LogFile     = log/${analysis}/${outFile}/${analysis}_${OutBase}.log
set ErrFile     = log/${analysis}/${outFile}/${analysis}_${OutBase}.err

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile

set arg = "$analysis $useEfficiency $triggerCoincidence $subLeadPtMin $leadPtMin $jetPtMax $jetRadius $constPtCut $binsEta $binsPhi $outLocation $outName $outNameTree $Files"

qsub -V -q erhiq -l mem=7GB -o $LogFile -e $ErrFile -N auauCorr -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

end
