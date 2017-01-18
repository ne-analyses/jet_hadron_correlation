# first make sure program is updated and exists
make bin/pythia_background || exit


set ExecPath = `pwd`
set analysis = $1
set execute = './bin/pythia_background'
set base = /nfs/rhi/STAR/Data/CleanAuAuY7/Clean

set nEvents = 1000

# Create the folder name for output
set outFile = pythia
# Make the directories since they may not exist...
if ( ! -d out/${outFile} ) then
mkdir -p out/${outFile}
endif

if ( ! -d log/pythia/${outFile} ) then
mkdir -p log/pythia/${outFile}
endif

# Now Submit jobs for each data file
foreach input ( ${base}* )

# Create the output file base name
set OutBase = `basename $input | sed 's/.root//g'`

# Make the output names and path
set outLocation = "out/${outFile}/"
set outName = corr_${OutBase}.root
set outNameTree = tree_${OutBase}.root

# Input files
set Files = ${input}

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out
set LogFile     = log/pythia/${outFile}/${OutBase}.log
set ErrFile     = log/pythia/${outFile}/${OutBase}.err

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile

set arg = " $nEvents $outLocation $outName $outNameTree $Files"

qsub -V -q erhiq -l mem=2GB -o $LogFile -e $ErrFile -N pythiaBkg -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

end
