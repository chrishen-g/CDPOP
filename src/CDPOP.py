# CDPOP.py
# Author: Erin L Landguth
# Created: February 2008
# v 1.2 Release: December 2011
# ----------------------------------------------------------------------------
# General CDPOP information
appName = "CDPOP"
appVers = "version 1.3.20"
appRele = "2024.09.12-09:55:00MDT"
authorNames = "Erin L Landguth et al."

# ---------------
# Global symbols
#----------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = True
# File absolute paths for importing functions
SRC_PATH =  "../src/"

# ------------------------------------------
# Import Modules with Except/Try statements
# ------------------------------------------
# Python specific functions
import datetime,time,pdb,os,sys,shutil
from tqdm.contrib.concurrent import process_map

import functools

# Numpy functions
try:
	import numpy as np                    
except ImportError as eMsg:
	print(("ImportError (%s) Numpy required."%(eMsg)))
	sys.exit(-1)

#Import the package specific folders
CDPOP_folder = os.path.dirname(os.path.abspath(SRC_PATH+"CDPOP"))

if CDPOP_folder not in sys.path:
     sys.path.insert(0, CDPOP_folder)

# CDPOP functions
try:
	from CDPOP_Modules import * 
except ImportError:
	raise ImportError("CDPOP_Modules required.")
try:
	from CDPOP_PostProcess import *
except ImportError:
	raise ImportError("CDPOP_PostProcess required.")
try:
	from CDPOP_PreProcess import *
except ImportError:
	raise ImportError("CDPOP_PreProcess required.")
try:
	from CDPOP_Mate import *
except ImportError:
	raise ImportError("CDPOP_Mate required.")
try:
	from CDPOP_Offspring import *
except ImportError:
	raise ImportError("CDPOP_Offspring required.")
try:
	from CDPOP_Disperse import *
except ImportError:
	raise ImportError("CDPOP_Disperse required.")
from CDPOP_montecarlo import mcrun

#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 
if __name__ == '__main__':
		
	# ------------------------------------------------------	
	# Start timer, get script arguments, create log writeout
	# ------------------------------------------------------
	# Timing events: start
	start_time = datetime.datetime.now()
	foldertime = int(time.time())
	
	if len(sys.argv) >= 4:
		datadir = sys.argv[1]+'/'
		fileans = datadir+sys.argv[2]
		outdir = datadir+sys.argv[3]+str(foldertime)+'/'
	
	# If user did not specify .rip file
	else:
		print("User must specify data directory, input file name, and output file directory (e.g., at command line type CDPOP.py ../CDPOP_data/ inputvariables16pnts.csv exampleout).")
		sys.exit(-1)	
	
	# If .ip file does not exist
	if not os.path.exists(fileans):
		print(("Cannot find or open runtime inputs file(%s)"%(fileans)))
		sys.exit(-1)
	
	# Create output file directory - will automatically put in the data directory
	os.mkdir(outdir)
	
	# This properly names log file
	logSessionPath = outdir+"cdpop.log"
	logfHndl =open(logSessionPath,'w')
	
	msgVerbose = True
	logMsg(logfHndl,"\n%s Release %s Version %s\n"%(appName,appRele,appVers))
	logMsg(logfHndl,"Author(s): %s"%(authorNames)+'\n')
	logMsg(logfHndl,"Session runtime inputs from: %s"%(fileans)+'\n\n')    
	msgVerbose = False
	
	# ------------------------------------	
	# Call DoUserInput()
	# ------------------------------------
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Call function and store inputvariables
	batchVars,batchVarsIndex,nSimulations = loadFile(fileans,1,',',True)
	
	# Print to log
	stringout = 'DoUserInput(): '+str(datetime.datetime.now() -start_time1) + ''
	logMsg(logfHndl,stringout)
	print('DoUserInput(): ',str(datetime.datetime.now() -start_time1),'')

	# -------------------------------------	
	# Begin Batch Looping
	# -------------------------------------
	# This loop is defined by the number of rows in inputvariables.csv
	for ibatch in range(nSimulations):
	
		# Timing events: start
		start_timeB = datetime.datetime.now()
		
		# Store all information and the type of each, also do some error checks 
		xyfilename = batchVars['xyfilename'][ibatch]
		allefreqfilename = batchVars['allefreqfilename'][ibatch]
		agefilename = batchVars['agefilename'][ibatch]
		matecdmatfile = batchVars['matecdmat'][ibatch]
		dispcdmatfile = batchVars['dispcdmat'][ibatch]
		mcruns = int(batchVars['mcruns'][ibatch])
		looptime = int(batchVars['looptime'][ibatch])
		nthfile_out = batchVars['output_years'][ibatch]
		cdclimgentimelist = batchVars['cdclimgentime'][ibatch]
		unicor_out = batchVars['output_unicor'][ibatch]		
		matemoveno = batchVars['matemoveno'][ibatch]
		matemoveparA = batchVars['matemoveparA'][ibatch]
		matemoveparB = batchVars['matemoveparB'][ibatch]
		matemoveparC = batchVars['matemoveparC'][ibatch]
		matemovethresh = batchVars['matemovethresh'][ibatch]
		freplace = batchVars['Freplace'][ibatch]
		mreplace = batchVars['Mreplace'][ibatch]
		mpaternity = batchVars['multiple_paternity'][ibatch]
		selfans = batchVars['selfans'][ibatch]
		matefreq = float(batchVars['mateFrequency'][ibatch])
		sexans = batchVars['sexans'][ibatch]
		Fdispmoveno = batchVars['Fdispmoveno'][ibatch]
		FdispmoveparA = batchVars['FdispmoveparA'][ibatch]
		FdispmoveparB = batchVars['FdispmoveparB'][ibatch]
		FdispmoveparC = batchVars['FdispmoveparC'][ibatch]
		Fdispmovethresh = batchVars['Fdispmovethresh'][ibatch]
		Mdispmoveno = batchVars['Mdispmoveno'][ibatch]
		MdispmoveparA = batchVars['MdispmoveparA'][ibatch]
		MdispmoveparB = batchVars['MdispmoveparB'][ibatch]
		MdispmoveparC = batchVars['MdispmoveparC'][ibatch]
		Mdispmovethresh = batchVars['Mdispmovethresh'][ibatch]
		philopatry = batchVars['philopatry'][ibatch]
		offnovals = batchVars['offno'][ibatch]
		Femalepercent = int(batchVars['Femalepercent'][ibatch])
		equalsexratio = batchVars['EqualsexratioBirth'][ibatch]
		twinning_pass = batchVars['TwinningPercent'][ibatch]
		popmodel = batchVars['popModel'][ibatch]
		K_envvals = batchVars['K_env'][ibatch]
		#subpopmort = batchVars['subpopmortperc'][ibatch]
		gridformat = batchVars['gridformat'][ibatch]
		muterate = float(batchVars['muterate'][ibatch])
		mutationans = batchVars['mutationtype'][ibatch]
		loci = int(batchVars['loci'][ibatch])
		intgenesans = batchVars['intgenesans'][ibatch]
		alleles = batchVars['alleles'][ibatch]
		mtdna = batchVars['mtdna'][ibatch]
		geneswap = int(batchVars['startGenes'][ibatch]) 
		cdevolveans = batchVars['cdevolveans'][ibatch]
		startSelection = int(batchVars['startSelection'][ibatch])
		betaFile_selection = batchVars['betaFile_selection'][ibatch]
		epigeneans = batchVars['epigeneans'][ibatch]
		startEpigene = int(batchVars['startEpigene'][ibatch])
		betaFile_epigene = batchVars['betaFile_epigene'][ibatch]
		cdinfect = batchVars['cdinfect'][ibatch]
		transmissionprob = float(batchVars['transmissionprob'][ibatch])
		matedist_out = batchVars['output_matedistance'][ibatch]
						
		# Distill and some error checking
		# -------------------------------
		# Grab the nthfile list range specific to user input, list or sequence
		if not isinstance(nthfile_out, (list,tuple)):
			nthfile_out = int(nthfile_out)
			if nthfile_out != 0:
				nthfile = list(range(0,looptime+nthfile_out,nthfile_out))
				del(nthfile[-1]) # Delete the last value 0, looptime - 1
			else:
				nthfile = [0]
		# If specified years with |
		else:
			nthfile = []
			# Split up list, removing space values, and appending to nthfile
			for inum in range(len(nthfile_out)):
				# Error check here if | at the end
				if len(nthfile_out[inum]) != 0:
					nthfile.append(int(nthfile_out[inum]))
		# Error check on nthfile, must be 1 less than looptime for indexing
		if max(nthfile) >= looptime:
			print('nthfile selection maximum value must be less than to looptime.')
			sys.exit(-1)
		
		# Store cdmat file information - header file (loadFile()) passes tuple or string if only 1
		if not isinstance(cdclimgentimelist, (list,tuple)):
			cdclimgentime = [cdclimgentimelist]
		else: 
			cdclimgentime = cdclimgentimelist
		if cdclimgentime[0] != '0':
			print('First cdclimate time must be 0.')
			sys.exit(-1)
	
		# Get mortality here: if tuple not returned and just one number applied across all adult ages 
		if not isinstance(agefilename,(list,tuple)):
			agefilename = [agefilename]
		else:
			agefilename = agefilename			
		if intgenesans == 'file' and allefreqfilename == 'N':
			print('Allele frequency file option specified, must give name of file.')
			sys.exit(-1)
		elif intgenesans == 'file_var' and allefreqfilename == 'N':
			print('Allele frequency file option specified, must give name of file.')
			sys.exit(-1)
		elif intgenesans == 'file' and allefreqfilename != 'N':
			if not isinstance(allefreqfilename,(list,tuple)):
				allefreqfilename = [allefreqfilename]
			else:
				allefreqfilename = allefreqfilename
		elif intgenesans == 'file_var' and allefreqfilename != 'N':
			if not isinstance(allefreqfilename,(list,tuple)):
				allefreqfilename = [allefreqfilename]
			else:
				allefreqfilename = allefreqfilename
		
		# If multiple XY files were specified for introducing individuals, then put in list as above
		if not isinstance(xyfilename,(list,tuple)):
			xyfilename = [xyfilename]
		else:	
			xyfilename = xyfilename
					
		# Create allele array
		if len(alleles.split(';')) == 1:
			alleles = int(batchVars['alleles'][ibatch])*np.ones(loci,int)
		else:
			alleles = np.asarray(alleles.split(';'),dtype = int)			
		
		# ---------------------------------
		# Some more Error checking
		# ---------------------------------
		# Have to have at least 2 alleles
		if len(np.where(alleles == 0)[0]) != 0:
			print('Must have at least 2 alleles per locus.')
			sys.exit(-1)
		if len(np.where(alleles == 1)[0]) != 0:
			print('Must have at least 2 alleles per locus.')
			sys.exit(-1)		
		
		# If cdevolve is turned on must have 2 alleles
		if cdevolveans != 'N' and alleles[0] != 2:
			print('Warning: More than 2 alleles per locus specified. CDEVOLVE only considers first 2 alleles in selection and epigenetic models, unless multiple loci models were specified.')
		if epigeneans != 'N' and alleles[0] != 2:
			print('Input Error: More than 2 alleles per locus specified. Epigenetics only considers 2 alleles in selection and epigenetic models.')
			sys.exit(-1)
		# Must have more than 1 loci
		if loci <= 1:
			print('Currently, CDPOP needs more than 1 locus to run.')
			sys.exit(-1)
			
		# Error check on forward mutation in A and backward mutation in B
		#	Can only happen if cdevolve == 2.
		if mutationans == 'forwardAbackwardBrandomN' and cdevolveans != '2':
			print('This special case of mutation is for AAbb ancestors and 2-locus selection.')
			sys.exit(-1)		
		
		# Check on parameters: equal sex ratio
		if (equalsexratio == 'N' or equalsexratio == 'AtBirth' or equalsexratio == 'WrightFisher') == False:
			print('Equal sex ratio parameter must be N, AtBirth, WrightFisher.')
			sys.exit(-1)
			
		# For female philopatry
		if philopatry == 'F' or philopatry == 'female' or philopatry == 'f' or philopatry == '0' or philopatry == 'Female': 
			
			if equalsexratio != 'AtBirth':
				print('Warning: Female philopatry is turned on and equal sex ratio at birth is recommended.')
			philopatry == 'F'
		# For male philopatry
		elif philopatry == 'M' or philopatry == 'male' or philopatry == 'm' or philopatry == '1' or philopatry == 'Male': 
			if equalsexratio != 'AtBirth':
				print('Warning: Male philopatry is turned on and equal sex ratio at birth is recommended.')
			philopatry == 'M'
		# Error if something else but N
		elif philopatry != 'N':
			print('Philopatry answer either has to be Male biased (M), Female biased (F), or unbiased (N).')
			sys.exit(-1)
					
		# grid format
		if (gridformat == 'cdpop' or gridformat == 'general' or gridformat == 'genalex' or gridformat == 'genepop' or gridformat == 'structure') == False:
			print('Grid format parameter not an option.')
			sys.exit(-1)
		
		# If genepop, some conditions
		if gridformat == 'genepop' and (len(np.where(alleles >= 99)[0]) > 0 or loci > 99):
			print('GENEPOP format requires less than 99 alleles and 99 loci.')
			sys.exit(-1)
			
		# For multiple paternity
		if mpaternity == 'Y' and (freplace != 'Y' or mreplace != 'Y'):
			print('Multiple paternity option is selected, then female and male with replacement must be both Y')
			sys.exit(-1)
		
		# Check burn in times
		if cdevolveans != 'N' and startSelection < geneswap:
			print('Start selection time must be less than genetic exchange start time (startGenes < startSelection).')
			sys.exit(-1)
		if epigeneans != 'N' and startEpigene < geneswap:
			print('Start epigenetics time must be less than genetic exchange start time (startGenes < startEpigene).')
			sys.exit(-1)
			
		# Check multiple selection model and number of loci
		if cdevolveans.split('_')[0] == 'M':
			if int(cdevolveans.split('_')[2].split('L')[1]) > loci:
				print('More loci under selection than specified number of total loci.')
				sys.exit(-1)
			if len(cdevolveans.split('_')) != 5:
				print('Multilocus selection specified, must have 6 arguments, see usermanual examples.')
				sys.exit(-1)
			if alleles[0] != int(cdevolveans.split('_')[3].split('A')[1]):
				print('Multilocus selection specified, must specify number of alleles for both neutral and selection markers.')
				sys.exit(-1)
			if intgenesans == 'random_var' or intgenesans == 'file_var':
				print('Variable allele assignment per locus can not be used with multi-locus selection.')
				sys.exit(-1)
		if epigeneans != 'N':
			if int(epigeneans.split('_')[1].split('L')[1]) > loci:
				print('More loci in epigenetic model than specified number of total loci.')
				sys.exit(-1)
			if len(epigeneans.split('_')) != 5:
				print('Epigenetic module specified, must have 4 arguments, see usermanual examples.')
				sys.exit(-1)
		if epigeneans != 'N' and cdevolveans.split('_')[0] == 'M':
			if int(cdevolveans.split('_')[2].split('L')[1]) + int(epigeneans.split('_')[1].split('L')[1]) > loci:
				print('More loci in epigenetic and selection model than specified number of total loci.')
				sys.exit(-1)		
		# Multiple files listed must equal cdclimgentime
		if len(xyfilename) > 1:
			if len(cdclimgentime) != len(xyfilename):
				print('Multiple xyfiles given, then must match the number of CDClimate generations.')
				sys.exit(-1)
			if intgenesans != 'file_introduce' and intgenesans != 'file_introduce_var':
				print('Multiple xyfiles given, then must use option file_introduce for genes.')
				sys.exit(-1)
			if len(allefreqfilename) != len(xyfilename):
				print('Multiple xyfiles given, then must match the number of allele frequency files given.')
				sys.exit(-1)
			if equalsexratio != 'N':
				print('Multiple xyfiles given, set equal sex ratio to N')
				sys.exit(-1)
			
		# Error checking for intgenesans special case for file_introduce
		if intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':
			if geneswap != 0:
				print('Gene swap starts immediately.')
		
		if popmodel == 'logistic':
			print('Logistic disabled temporarily, email erin.landguth@mso.umt.edu.')
			sys.exit(-1)
			
		# For Hindex answer
		if cdevolveans.split('_')[0] == 'Hindex':
			# Split for Gaussian
			if cdevolveans.split('_')[1] == 'Linear':
				if len(cdevolveans.split('_')[2].split(';')) != 6:
					print('CDEVOLVE answer is Hindex and 6 parameters for the Linear function must be specified, see user manual and example files.')
					sys.exit(-1)
			else:
				print('CDEVOLVE answer Hindex and only Linear option is allowed for now.')
				sys.exit(-1)
		
		# ---------------------------------------------	
		# Begin Monte-Carlo Looping
		# ---------------------------------------------
		
		# xrange(mcruns) is typically 10 - 50...and it takes a long time.
		# for ithmcrun in range(mcruns):

		args = (
			xyfilename, allefreqfilename, agefilename, matecdmatfile, dispcdmatfile,
			mcruns, looptime, nthfile_out, cdclimgentimelist, unicor_out,
			matemoveno, matemoveparA, matemoveparB, matemoveparC, matemovethresh,
			freplace, mreplace, mpaternity, selfans, matefreq, sexans,
			Fdispmoveno, FdispmoveparA, FdispmoveparB, FdispmoveparC, Fdispmovethresh,
			Mdispmoveno, MdispmoveparA, MdispmoveparB, MdispmoveparC, Mdispmovethresh,
			philopatry, offnovals, Femalepercent, equalsexratio, twinning_pass,
			popmodel, K_envvals, gridformat, muterate, mutationans, loci,
			intgenesans, alleles, mtdna, geneswap, cdevolveans, startSelection,
			betaFile_selection, epigeneans, startEpigene, betaFile_epigene,
			cdinfect, transmissionprob, matedist_out, datadir, fileans, outdir, ibatch, nthfile, cdclimgentime
		)

		mcrun_init = functools.partial(mcrun, args)

		process_map(mcrun_init, list(range(mcruns)))
		# End::Monte Carlo Loop
		
		# Print to log
		stringout = 'End Batch Loop'+str(ibatch)+': '+str(datetime.datetime.now() -start_timeB) + '\n'
		logMsg(logfHndl,stringout)
		print('End Batch Loop',str(ibatch),': ',str(datetime.datetime.now() -start_timeB),'\n')
		
	#End::Batch Loop
	
# End::Main Loop	
# Print to log
stringout = 'Total CDPOP Simulation Time: '+str(datetime.datetime.now() -start_time) + ''
logMsg(logfHndl,stringout)
logfHndl.close()
print('Total CDPOP Simulation Time: ',str(datetime.datetime.now() -start_time),'')