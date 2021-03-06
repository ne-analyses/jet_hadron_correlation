 #!/bin/csh

#  A Slightly ugly way to make it as simple as possible to run over different years,
#  And to make it easier to change variables and keep the book keeping simple
#  Command Line Arugments:
#  [1]: the year. Either y7 or y11
#  [2]: if set to 'default', uses original A_j parameterization
#  [2]: if not using default settings, sets whether a trigger is required to
#       be coincident with the leading jet
#
#  -----Di_jet selection criteria-----
#  [3]: leading jet minimum pt         - these set the minimum required  pt for 
#  [4]: subleading jet minimum pt     - hard dijet pairs
#  [5]: jet radius: set the radius in the jet definition (always with anti-kt algorithm)
#  [6]: hard constituent cut: the minimum pt for a particle to be clustered to find
#       hard di-jets
#  [7]: soft constituent cut: the minimum pt for a particle to be clustered to find
#       jets when we include the soft particles
#  [8]: number of bins to use in Vertex Z position
#  [9]: range to accept in Vz ( it is split symmetrically about Vz = 0 )
#  NOTE: [9]/[8] must be an integer or else...
#  [10]:whether or not to apply an efficiency correction track-by-track 
#
#  Output names and locations are generated by the script, and correspond to the above
#  variables

set year = $1
if ( $year == 'y11' ) then
    # Make sure the program exists
    make bin/dijet_hadron_correlation_y11 || exit
    # The data file base name
    set base = /Volumes/Promise\ Pegasus/STAR/Data/HaddedAuAu11picoNPE25_150526/AuAu11
    # Hot tower list used by reader
    set hotList = /Users/nickelsey/physics/software/event_structuredAu/badTowerList_y11.txt
    # Setting the command that will be executed
    set command = "./bin/dijet_hadron_correlation_y11"
else if ( $year == 'y7' ) then
    make bin/dijet_hadron_correlation_y7 || exit
    set base = /Users/nickelsey/physics/Data/SmallAuAu/Small
    set hotList = /Users/nickelsey/physics/software/eventStructuredAu/y7_AuAu_HT_hot_list.txt
    set command = "./bin/dijet_hadron_correlation_y7"
else
    echo 'Error: need to choose data set'
    exit
endif

if ( $# != "10" && !( $2 == 'default' ) ) then
    echo 'Error: illegal number of parameters'
    exit
endif

# Start the Condor File
echo "" > CondorFile
echo "Universe    = vanilla" >> CondorFile
echo "Executable  = ${command}" >> CondorFile
echo "getenv = true" >> CondorFile

# Arguments
set triggerCoincidence = $2
set leadPtMin = $3
set subLeadPtMin = $4
set jetRadius = $5
set hardConstituent = $6
set softConstituent = $7
set nVzBins = $8
set VzRange = $9
set useEfficiency = $10
    
if ( $2 == 'default' ) then
    set triggerCoincidence = 0
    set leadPtMin = 20.0
    set subLeadPtMin = 10.0
    set jetRadius = 0.4
    set hardConstituent = 2.0
    set softConstituent = 0.2
    set nVzBins = 12
    set VzRange = 60
    set useEfficiency = 0
endif

# Create the folder name for output
set outFile = lead_${leadPtMin}_sublead_${subLeadPtMin}_hard_${hardConstituent}_soft_${softConstituent}_vz_${nVzBins}_${VzRange}_useEff_${useEfficiency}
# Make the directories since they may not exist...
if ( ! -d out/${year}/${outFile} ) then
    mkdir -p out/${year}/${outFile}
    mkdir -p out/${year}/${outFile}/correlations
    mkdir -p out/${year}/${outFile}/dijets
    mkdir -p out/${year}/${outFile}/mixing
endif

# Same for logs
if ( ! -d logs/${year}/${outFile} ) then
    mkdir -p logs/${year}/${outFile}
endif

# Write the information for mixing to file
# Will be written to  the output location
# for the event mixing
set mixPath = out/${year}/${outFile}/mixing
set mixFile = ${mixPath}/mixing_params
echo $nVzBins > $mixFile
echo $VzRange >> $mixFile
echo $useEfficiency >> $mixFile

# Now Submit jobs for each data file
foreach input ( ${base}* )

    # Create the output file base name
    set OutBase = `basename $input | sed 's/.root//g'`

    # Make the output names and path
    set outLocation = "out/${year}/${outFile}/"
    set outName = correlations/corr_${OutBase}.root
    set outNameDiJets = dijets/dijets_${OutBase}.root

    # Input files
    set Files = ${input}

    # Logfiles. Thanks cshell for this "elegant" syntax to split err and out
    set LogFile     = logs/${year}/${outFile}/${year}_${OutBase}.out
    set ErrFile     = logs/${year}/${outFile}/${year}_${OutBase}.err

    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile

    set arg = "$triggerCoincidence $leadPtMin $subLeadPtMin $jetRadius $hardConstituent $softConstituent $nVzBins  $VzRange $useEfficiency $outLocation $outName $outNameDiJets $hotList $Files"
    set execute = "${command} ${arg}"

    # Write to CondorFile
    echo "Executing " $execute

    echo " " >> CondorFile
    echo "Output    = ${LogFile}" >> CondorFile
    echo "Error     = ${ErrFile}" >> CondorFile
    echo "Arguments = ${arg}" >> CondorFile
    echo "Queue" >> CondorFile

    
end

# Submit everything to condor
#condor_submit CondorFile






