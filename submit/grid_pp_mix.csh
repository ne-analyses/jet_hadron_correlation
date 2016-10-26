#!/bin/csh

# Script for event mixing auau trees on the grid
# Command line arguments:
# [1]: input directory ( it should have a long string of analysis settings generated from grid_auau_corr.csh )
# [2]: event mixing data ( either .root file or list )
# [3]: Is the data MB or HT?
# [4]: total number of events to look through
# [5]: number of events to mix with each trigger
#
# Can set default settings by only giving
# [1]: input directory
# [2]: 'default'


# first make sure program is updated and exists
make bin/event_mixing || exit

set ExecPath = `pwd`
set inputDir = $1
set execute = './bin/event_mixing'
set base = ${inputDir}/tree/tree

if ( $# != "5" && !( $2 == 'default' ) ) then
echo 'Error: illegal number of parameters'
exit
endif


# define arguments
set mixEvents = $2
set dataType = $3
set nEvents = $4
set eventsPerTrigger = $5

# make a base directory for logging
set logBase = `basename $inputDir`

#made the log directory
if ( ! -d log/mix/${logBase} ) then
mkdir -p log/mix/${logBase}
endif

if ( $2 == 'default' ) then

set mixEvents = 'pp_list/grid/evPP.list'
set dataType = 'HT'
set nEvents = '-1'
set eventsPerTrigger = '4000'
endif

# Now Submit jobs for each data file
foreach input ( ${base}*.root )

# Create the output file base name
set OutBase = `basename $input | sed 's/.root//g'`

# Make the output names and path
set outName = mixing/mix_${OutBase}.root

# Input files
set Files = ${input}

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out
set LogFile     = log/mix/${logBase}/mix_${OutBase}.log
set ErrFile     = log/mix/${logBase}/mix_${OutBase}.err

# get relative tree location
set treeFile = `basename $input`
set relativeTreeFile = tree/${treeFile}

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile

set arg = "$inputDir $relativeTreeFile $outName $dataType $nEvents $eventsPerTrigger $mixEvents"

qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ppMix -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

end
