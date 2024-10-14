import datetime
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


def mcrun(args, ithmcrun):
	# Timing events: start
	start_timeMC = datetime.datetime.now()

	(
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
	) = args

	# -----------------------------------------
	# Create storage variables
	# ------------------------------------------
	# These variables will be stored in output.csv at the end of the simulation
	Population = []
	Population_age = []
	Females_age = []
	Males_age = []
	Migrants = []
	Open = []
	Track_MDeaths = []
	Track_FDeaths = []
	Births = []
	Track_MOffDeaths = []
	Track_FOffDeaths = []
	DisperseDeaths = []
	CouldNotDisperse = []
	Opt3SelectionDeaths = []
	ToTFemales = []
	ToTMales = []
	BreedFemales = []
	BreedFemales_age = []
	BreedMales = []
	Female_BreedEvents = []
	females_nomate = []
	males_nomate = []
	Alleles = []
	He = []
	Ho = []
	AllelesMutated = []
	MateDistED = []
	FDispDistED = []
	MDispDistED = []
	MateDistCD = []
	FDispDistCD = []
	MDispDistCD = []
	MateDistEDstd = []
	FDispDistEDstd = []
	MDispDistEDstd = []
	MateDistCDstd = []
	FDispDistCDstd = []
	MDispDistCDstd = []
	Infected = []
	p1 = []
	p2 = []
	q1 = []
	q2 = []
	subpopmigration = []
	subpopemigration = []
	FAvgMate = []
	MAvgMate = []
	FSDMate = []
	MSDMate = []
	MateDistances = []
	Twins = []
	Track_EpigeneMod1 = []
	Track_EpigeneMod2 = []
	Track_EpigeneDeaths = []
	Track_EpigeneReset1 = []
	Track_EpigeneReset2 = []
	maxfit = []
	minfit = []

	# ------------------------------------
	# Call DoPreProcess()
	# ------------------------------------

	# Timing events: start
	start_time1 = datetime.datetime.now()

	# Call function
	tupPreProcess = DoPreProcess(outdir, ibatch, ithmcrun, \
								 xyfilename, agefilename, equalsexratio, loci, intgenesans, allefreqfilename,
								 alleles, 0, cdevolveans, cdinfect, Infected, \
								 subpopmigration, subpopemigration, datadir, geneswap, epigeneans, unicor_out)

	ithmcrundir = tupPreProcess[0]
	FID = tupPreProcess[1]
	id = tupPreProcess[2]
	sex = tupPreProcess[3]
	age = tupPreProcess[4]
	xgrid = tupPreProcess[5]
	xgridcopy = copy.deepcopy(xgrid)
	ygrid = tupPreProcess[6]
	ygridcopy = copy.deepcopy(ygrid)
	genes = tupPreProcess[7]
	nogrids = tupPreProcess[8]
	subpop = tupPreProcess[9]
	fitvals_pass = tupPreProcess[10]
	infection = tupPreProcess[11]
	Infected = tupPreProcess[12]
	subpopmigration = tupPreProcess[13]
	subpopemigration = tupPreProcess[14]
	Magemortvals = tupPreProcess[15]
	Fagemortvals = tupPreProcess[16]
	egg_lmbdavals = tupPreProcess[17]
	egg_sigmavals = tupPreProcess[18]
	allelst = tupPreProcess[19]
	Mnewmortperc = tupPreProcess[20]
	Fnewmortperc = tupPreProcess[21]
	Mmaturevals = tupPreProcess[22]
	Fmaturevals = tupPreProcess[23]
	intgenesans = tupPreProcess[24]
	xvars_betas_pass = tupPreProcess[25]
	epimod_pass = tupPreProcess[26]
	epireset_pass = tupPreProcess[27]
	hindex = tupPreProcess[28]
	gridmort_pass = tupPreProcess[29]

	# Print to log
	stringout = 'DoPreProcess(): ' + str(datetime.datetime.now() - start_time1) + ''
	# logMsg(logfHndl, stringout)
	print('DoPreProcess(): ', str(datetime.datetime.now() - start_time1), '')

	# -------------------------------------------
	# Start Generation Looping
	# -------------------------------------------

	# Begin generation loop
	for gen in range(looptime):

		# Timing events: start
		start_timeGen = datetime.datetime.now()

		# ---------------------------------
		# Call CDClimate()
		# ---------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()

		# Check gen time equal to cdclimgentime
		if len(np.where(np.asarray(cdclimgentime) == str(gen))[0]) == 1:
			tupClimate = DoCDClimate(datadir, np.where(np.asarray(cdclimgentime) == str(gen))[0][0], cdclimgentime,
									 matecdmatfile, dispcdmatfile, matemoveno, Fdispmoveno, Mdispmoveno,
									 matemovethresh, Fdispmovethresh, Mdispmovethresh, matemoveparA, matemoveparB,
									 matemoveparC, FdispmoveparA, FdispmoveparB, FdispmoveparC, MdispmoveparA,
									 MdispmoveparB, MdispmoveparC, subpop, Magemortvals, Fagemortvals, offnovals,
									 egg_lmbdavals, egg_sigmavals, K_envvals, Mnewmortperc, Fnewmortperc,
									 fitvals_pass, twinning_pass, Mmaturevals, Fmaturevals, betaFile_selection,
									 xvars_betas_pass, epimod_pass, epireset_pass, betaFile_epigene, cdevolveans,
									 epigeneans, gridmort_pass)

			cdmatrix_mate = tupClimate[0]
			cdmatrix_F = tupClimate[1]
			cdmatrix_M = tupClimate[2]
			thresh_mate = tupClimate[3]
			thresh_F = tupClimate[4]
			thresh_M = tupClimate[5]
			Fdisp_ScaleMin = tupClimate[6]
			Fdisp_ScaleMax = tupClimate[7]
			Mdisp_ScaleMin = tupClimate[8]
			Mdisp_ScaleMax = tupClimate[9]
			mate_ScaleMin = tupClimate[10]
			mate_ScaleMax = tupClimate[11]
			Magemort = tupClimate[12]
			Fagemort = tupClimate[13]
			offno = tupClimate[14]
			eggs_lambda = tupClimate[15]
			eggs_sigma = tupClimate[16]
			K_env = tupClimate[17]
			Mnewmort = tupClimate[18]
			Fnewmort = tupClimate[19]
			fitvals = tupClimate[20]
			mateno = tupClimate[21]
			Fdispno = tupClimate[22]
			Mdispno = tupClimate[23]
			twinning = tupClimate[24]
			Mmature = tupClimate[25]
			Fmature = tupClimate[26]
			betas_selection = tupClimate[27]
			xvars_betas = tupClimate[28]
			epimod = tupClimate[29]
			epireset = tupClimate[30]
			betas_epigene = tupClimate[31]
			gridmort = tupClimate[32]

			# Error check for if nofiles == nogrids system exit
			if nogrids != len(cdmatrix_mate):
				print('The cost distance matrix dimensions are not the same as the number of individuals.')
				sys.exit(-1)

			# Print to log
			stringout = 'DoCDCliamte(): ' + str(datetime.datetime.now() - start_time1) + ''
			# logMsg(logfHndl, stringout)
			print('DoCDClimate(): ', str(datetime.datetime.now() - start_time1), '')

		# -------------------------------
		# Call ReadGrid()
		# -------------------------------
		# Use information generated from PreProcess step for first
		#	generation, else use the following updated grid information
		if gen != 0:

			# Timing events: start
			start_time1 = datetime.datetime.now()
			'''move to before disperse
			# ---------------------------------
			# Add individuals if specified
			# ---------------------------------					
			if len(xyfilename) > 1:
				# Check gen time equal to cdclimgentime
				if len(np.where(np.asarray(cdclimgentime) == str(gen))[0]) == 1:

					# Timing events: start
					start_time1 = datetime.datetime.now()

					tupAddInds = AddIndividuals(cdclimgentime,gen,idnew,agenew,genesnew,sexnew,subpopnew,infectionnew,allelst,xyfilename,datadir,alleles,hindexnew)

					idnew = tupAddInds[0]
					sexnew = tupAddInds[1]
					agenew = tupAddInds[2]
					genesnew = tupAddInds[3]
					infectionnew = tupAddInds[4]
					hindexnew = tupAddInds[5]

					# Print to log
					stringout = 'AddIndividuals(): '+str(datetime.datetime.now() -start_time1) + ''
					logMsg(logfHndl,stringout)
					print(('AddIndividuals()',str(datetime.datetime.now() -start_time1),''))		
			'''

			tupReadGrid = ReadGrid(FIDnew, idnew, agenew, xgridnew, \
								   ygridnew, genesnew, equalsexratio, sexnew, subpopnew, \
								   infectionnew, allelst, geneswap, gen, intgenesans, hindexnew)

			FID = tupReadGrid[0]
			sex = tupReadGrid[1]
			id = tupReadGrid[2]
			age = tupReadGrid[3]
			xgrid = tupReadGrid[4]
			xgridcopy = tupReadGrid[5]
			ygrid = tupReadGrid[6]
			ygridcopy = tupReadGrid[7]
			genes = tupReadGrid[8]
			nogrids = tupReadGrid[9]
			subpop = tupReadGrid[10]
			infection = tupReadGrid[11]
			filledgrids = tupReadGrid[12]
			hindex = tupReadGrid[13]

			# Exit system if population is 0 or 1
			if filledgrids == 0 or filledgrids == 1:
				stringout = 'Population went extinct, program ended.'
				# logMsg(logfHndl, stringout)
				print(('Population went extinct after generation ' + str(gen - 1) + '.\n'))
				break
			# Here we system exit if there are only F or M left
			if sexans == 'Y' and len(np.where(np.asarray(sex) == '0')[0]) == 0:
				# Print to log
				stringout = 'No females left in population, program ended.'
				# logMsg(logfHndl, stringout)
				print(('There are no more females left in population after generation ' + str(gen - 1) + '.\n'))
				break
			if sexans == 'Y' and len(np.where(np.asarray(sex) == '1')[0]) == 0:
				# Print to log
				stringout = 'No males left in population, program ended.'
				# logMsg(logfHndl, stringout)
				print(('There are no more males left in population after generation ' + str(gen - 1) + '.\n'))
				break
			# Print to log
			stringout = 'ReadGrid(): ' + str(datetime.datetime.now() - start_time1) + ''
			# logMsg(logfHndl, stringout)
			print('ReadGrid(): ', str(datetime.datetime.now() - start_time1), '')

		# ---------------------------------
		# Call GetMetrics()
		# ---------------------------------
		# pdb.set_trace()
		# Timing events: start
		start_time1 = datetime.datetime.now()

		tupGetMetrics = GetMetrics(Population, nogrids, loci, alleles, genes, \
								   gen, Ho, Alleles, He, subpop, p1, p2, q1, q2, Population_age, Females_age,
								   Males_age, age, sex, Magemort, geneswap, cdevolveans, xvars_betas,
								   betas_selection, maxfit, minfit)

		filledgrids = tupGetMetrics[0]
		subgridtotal = tupGetMetrics[1]

		# Print to log
		stringout = 'GetMetrics(): ' + str(datetime.datetime.now() - start_time1) + ''
		# logMsg(logfHndl, stringout)
		print('GetMetrics(): ', str(datetime.datetime.now() - start_time1), '')
		# pdb.set_trace()
		# --------------------------------
		# Call DoEpigenetics
		# --------------------------------
		# Add time step to each tracker variable here
		if epigeneans != 'N':
			# Timing events: start
			start_time1 = datetime.datetime.now()

			tupEpigene = DoEpigenetics(epimod, betas_epigene, sex, id, age, genes, infection, Track_EpigeneMod1,
									   Track_EpigeneMod2, Track_EpigeneDeaths, gen, cdevolveans, epigeneans,
									   startEpigene, geneswap, alleles, loci)

			id = tupEpigene[0]
			sex = tupEpigene[1]
			age = tupEpigene[2]
			genes = tupEpigene[3]
			infection = tupEpigene[4]

			# Print to log
			stringout = 'DoEpigenetics(): ' + str(datetime.datetime.now() - start_time1) + ''
			# logMsg(logfHndl, stringout)
			print('DoEpigenetics(): ', str(datetime.datetime.now() - start_time1), '')

		# ---------------------------------------
		# Call DoMate()
		# ---------------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()

		tupMate = DoMate(nogrids, sex, age, \
						 freplace, mreplace, mateno, thresh_mate, \
						 cdmatrix_mate, MateDistED, MateDistCD, xgridcopy, \
						 ygridcopy, ToTMales, ToTFemales, BreedMales, BreedFemales, \
						 sexans, selfans, \
						 MateDistEDstd, MateDistCDstd, FAvgMate, MAvgMate, \
						 FSDMate, MSDMate, filledgrids, Female_BreedEvents, gen, subpop, BreedFemales_age, Magemort,
						 Mmature, Fmature, mate_ScaleMax, mate_ScaleMin, matemoveparA, matemoveparB, matemoveparC,
						 MateDistances, matefreq)
		Bearpairs = tupMate[0]
		females = tupMate[1]
		females_nomate.append(tupMate[2])
		males = tupMate[3]
		males_nomate.append(tupMate[4])
		mature = tupMate[5]

		# Temporary, get rid of this eventually...
		CDpairs = copy.deepcopy(Bearpairs)

		# Print to log
		stringout = 'DoMate(): ' + str(datetime.datetime.now() - start_time1) + ''
		# logMsg(logfHndl, stringout)
		print('DoMate(): ', str(datetime.datetime.now() - start_time1), '')

		# ---------------------------------------
		# Call DoOffspring()
		# ---------------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()

		tupDoOff = DoOffspring(offno, eggs_lambda, Bearpairs, CDpairs, Femalepercent, \
							   Births, infection, transmissionprob, equalsexratio, \
							   Mnewmort, Fnewmort, Track_MOffDeaths, Track_FOffDeaths, eggs_sigma, age, sex,
							   twinning, Twins, subpop)

		offspring = tupDoOff[0]
		offspringno = tupDoOff[1]

		# Print to log
		stringout = 'DoOffspring(): ' + str(datetime.datetime.now() - start_time1) + ''
		# logMsg(logfHndl, stringout)
		print('DoOffspring(): ', str(datetime.datetime.now() - start_time1), '')

		# ---------------------------------------
		# Call InheritGenes()
		# ---------------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()

		offspring = InheritGenes(gen, AllelesMutated, offspringno, offspring, genes, loci, muterate, mtdna, \
								 mutationans, geneswap, epireset, Track_EpigeneReset1, Track_EpigeneReset2, \
								 startEpigene, epigeneans, cdevolveans, alleles, hindex)

		# Print to log
		stringout = 'InheritGenes(): ' + str(datetime.datetime.now() - start_time1) + ''
		# logMsg(logfHndl, stringout)
		print('InheritGenes(): ', str(datetime.datetime.now() - start_time1), '')

		# ------------------------------------------
		# Call DoAdultMortality()
		# ------------------------------------------
		# pdb.set_trace()
		# Timing events: start
		start_time1 = datetime.datetime.now()

		tupAMort = DoMortality(nogrids, sex, id, age, xgrid, ygrid, gen, genes, Track_MDeaths, Track_FDeaths, \
							   FID, Magemort, Fagemort, infection, geneswap, popmodel, K_env, fitvals, \
							   mature, cdevolveans, Opt3SelectionDeaths, startSelection, subpop, hindex)

		freegrid = tupAMort[0]
		id = tupAMort[1]
		sex = tupAMort[2]
		age = tupAMort[3]
		xgrid = tupAMort[4]
		ygrid = tupAMort[5]
		genes = tupAMort[6]
		FID = tupAMort[7]
		infection = tupAMort[8]  # coming out str, change to int
		# Mature skipped, not used from here on out
		hindex = tupAMort[10]

		# Print to log
		stringout = 'DoAdultMortality(): ' + str(datetime.datetime.now() - start_time1) + ''
		# logMsg(logfHndl, stringout)
		print('DoAdultMortality(): ', str(datetime.datetime.now() - start_time1), '')
		# pdb.set_trace()
		# ----------------------------------------------
		# Add individuals if specified, starting gen = 1
		# ----------------------------------------------
		if len(xyfilename) > 1 and gen != 0:
			# Check gen time equal to cdclimgentime
			if len(np.where(np.asarray(cdclimgentime) == str(gen))[0]) == 1:
				# Timing events: start
				start_time1 = datetime.datetime.now()

				# tupAddInds = AddIndividuals(cdclimgentime,gen,id,age,genes,sexnew,subpop,infection,allelst,xyfilename,datadir,alleles,hindex)
				tupAddInds = AddIndividuals(cdclimgentime, gen, id.tolist(), age.tolist(), genes, sex.tolist(),
											FID.tolist(), infection.tolist(), allelst, xyfilename, datadir, alleles,
											hindex.tolist(), subpop, freegrid)

				id = tupAddInds[0]
				sex = tupAddInds[1]
				age = tupAddInds[2]
				genes = tupAddInds[3]
				infection = tupAddInds[4]
				hindex = tupAddInds[5]
				FID = tupAddInds[6]
				freegrid = tupAddInds[7]

				# Print to log
				stringout = 'AddIndividuals(): ' + str(datetime.datetime.now() - start_time1) + ''
				# logMsg(logfHndl, stringout)
				print(('AddIndividuals()', str(datetime.datetime.now() - start_time1), ''))
		# pdb.set_trace()
		# ------------------------------------------
		# Call DoDisperse()
		# ------------------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()

		tupDoDisp = DoDisperse(offspringno, freegrid, offspring, Fdispno, \
							   Mdispno, cdmatrix_F, cdmatrix_M, gen, \
							   Migrants, Open, loci, alleles, \
							   xgridcopy, ygridcopy, FDispDistED, MDispDistED, FDispDistCD, MDispDistCD, \
							   cdevolveans, fitvals, FDispDistEDstd, MDispDistEDstd, \
							   FDispDistCDstd, MDispDistCDstd, subpop, subpopmigration, DisperseDeaths,
							   CouldNotDisperse, \
							   gridmort, philopatry, females, subpopemigration, females_nomate[gen], \
							   males, males_nomate[gen], startSelection, thresh_F, thresh_M, Fdisp_ScaleMax, \
							   Fdisp_ScaleMin, Mdisp_ScaleMax, Mdisp_ScaleMin, FdispmoveparA, FdispmoveparB, \
							   FdispmoveparC, MdispmoveparA, MdispmoveparB, MdispmoveparC, betas_selection,
							   xvars_betas, maxfit, minfit)

		OffDisperseIN = tupDoDisp[0]
		opengrids = tupDoDisp[1]
		DispDistCD = tupDoDisp[2]

		# Print to log
		stringout = 'DoDisperse(): ' + str(datetime.datetime.now() - start_time1) + ''
		# logMsg(logfHndl, stringout)
		print('DoDisperse(): ', str(datetime.datetime.now() - start_time1), '')

		# ------------------------------------------
		# Call DoOutput()
		# ------------------------------------------
		# pdb.set_trace()
		# Timing events: start
		start_time1 = datetime.datetime.now()

		tupDoOut = DoOutput(nogrids, FID, OffDisperseIN, \
							xgridcopy, ygridcopy, gen, id, sex, age, xgrid, \
							ygrid, genes, nthfile, ithmcrundir, loci, alleles, subpop, \
							gridformat, infection, Infected, cdinfect, \
							opengrids, DispDistCD, geneswap, hindex, unicor_out)

		FIDnew = tupDoOut[0]
		idnew = tupDoOut[1]
		sexnew = tupDoOut[2]
		agenew = tupDoOut[3]
		xgridnew = tupDoOut[4]
		ygridnew = tupDoOut[5]
		genesnew = tupDoOut[6]
		subpopnew = tupDoOut[7]
		infectionnew = tupDoOut[8]  # mix of string/ints, fix
		hindexnew = tupDoOut[9]

		# Print to log
		stringout = 'DoOutput(): ' + str(datetime.datetime.now() - start_time1) + ''
		# logMsg(logfHndl, stringout)
		print('DoOutput(): ', str(datetime.datetime.now() - start_time1), '')

		# Print to log
		stringout = 'End Generation Loop' + str(gen) + ': ' + str(datetime.datetime.now() - start_timeGen) + '\n'
		# logMsg(logfHndl, stringout)
		print('End Generation Loop', str(gen), ': ', str(datetime.datetime.now() - start_timeGen), '\n')

	# End::generation loop

	# ------------------------------------------
	# Call DoPostProcess()
	# ------------------------------------------

	# Timing events: start
	start_time1 = datetime.datetime.now()

	DoPostProcess(ithmcrundir, nogrids, \
				  xgridcopy, ygridcopy, \
				  loci, alleles, looptime, Population, ToTFemales, ToTMales, \
				  BreedFemales, BreedMales, Migrants, Births, \
				  Track_MDeaths, Track_FDeaths, Alleles, He, Ho, AllelesMutated, \
				  MateDistED, FDispDistED, MDispDistED, MateDistCD, FDispDistCD, MDispDistCD, nthfile, \
				  p1, p2, q1, q2, Infected, subpop, MateDistEDstd, \
				  FDispDistEDstd, MDispDistEDstd, MateDistCDstd, FDispDistCDstd, MDispDistCDstd, subpopmigration, \
				  FAvgMate, MAvgMate, FSDMate, MSDMate, DisperseDeaths, Open, CouldNotDisperse, \
				  Female_BreedEvents, gridformat, subpopemigration, females_nomate, subgridtotal, Track_MOffDeaths,
				  Track_FOffDeaths, Population_age, Females_age, Males_age, \
				  BreedFemales_age, Opt3SelectionDeaths, MateDistances, matedist_out, Twins, Track_EpigeneMod1,
				  Track_EpigeneMod2, Track_EpigeneDeaths, Track_EpigeneReset1, Track_EpigeneReset2)

	# Print to log
	stringout = 'DoPostProcess(): ' + str(datetime.datetime.now() - start_time1) + ''
	# logMsg(logfHndl, stringout)
	print('DoPostProcess(): ', str(datetime.datetime.now() - start_time1), '')

	# Print to log
	stringout = 'End Monte Carlo Loop' + str(ithmcrun) + ': ' + str(datetime.datetime.now() - start_timeMC) + '\n'
	# logMsg(logfHndl, stringout)
	print('End Monte Carlo Loop', str(ithmcrun), ': ', str(datetime.datetime.now() - start_timeMC), '\n')