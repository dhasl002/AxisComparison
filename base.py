import chimera
from chimera.baseDialog import ModelessDialog
from chimera import help, openModels, Molecule
from chimera import viewer
from Group import Group, GroupAttr
import Tix, Pmw
import Tkinter
import CGLtk
import os
import chimera
import subprocess
from chimera import dialogs
from Tkinter import Tk
import Queue
import multiprocessing
import ModelPanel
import platform
import time
import shutil
from VolumeViewer import Volume

_buttonInfo = {}
_mp = None
oPath = ""

def addButton(name, callback, minModels=1, maxModels=None,
			moleculesOnly=True, balloon=None, defaultFrequent=True,
			defaultFavorite=None, groupsOkay=False):
	"""Add a button to the 'Model Actions' button list.

	   'name' is the button name (duh).  'callback' is the
	   function to call when the button is pressed.  The arg to
	   'callback' will be a list of models.  'min/maxModels'
	   indicate how many models have to be selected in the
	   browser for the button to be active ('None' indicates no
	   limit).  if 'moleculesOnly' is True, then those models have
	   to be Molecules.

	   This is a module function so that it can be called even if
	   the model panel has not yet been created.
	"""

	if defaultFavorite is None:
		defaultFavorite = defaultFrequent
	if _buttonInfo.has_key(name):
		raise KeyError, \
			"Button named '%s' already exists" % name
	_buttonInfo[name] = (callback, minModels, maxModels,
					moleculesOnly, balloon, defaultFavorite, groupsOkay)

	if _mp:
		_mp._confDialog.newButton(name, balloon=balloon,
					defaultFavorite=defaultFavorite)
		_mp._showButton(name)

_columnNames = []
_valueTypes = []
_valueFuncs = []
_defaultShows = []
def addColumns(columnInfo, defaultShown=1):
	"""Add columns to the model table.

	   'columnInfo' is a list of 3-tuples, one for each column
	   to add.  The tuple consists of (column name, value type,
	   value-fetch function).  The value type should be 'text',
	   'image', 'imagetext', or 'toggle'.  The value-fetch function
	   takes one argument (a model) and (for 'image' and 'text')
	   should return the value to display in the table cell.  For
	   'imagetext' the return value should be an (image, text)
	   tuple.  'toggle' shows a toggle button and the return value
	   should be a (boolean, callback) tuple.  If the boolean is
	   true, a check will be shown on the toggle button; otherwise
	   the button is blank.  The callback is invoked when the
	   toggle is pressed, with the model and the new boolean as
	   arguments.  The value of an image is the name of the image
	   (the Tix name, e.g. 'tick' for tickmark).  A value of None
	   for image or text will leave a blank cell.

	   'defaultShown' controls whether the column is shown in
	   the model table or not as long as the user has not yet
	   expressed a preference in the Configuration panel about it.
	"""

	noneShown = 1
	for name,type,func in columnInfo:
		if name in _columnNames:
			raise ValueError, "Duplicate model panel"\
						"column name: %s" % name
		_columnNames.append(name)
		_valueTypes.append(type)
		_valueFuncs.append(func)
		_defaultShows.append(defaultShown)
		if _mp:
			try:
				shown = _mp._confDialog.prefs[
							'shownColumns'][name]
			except KeyError:
				shown = defaultShown
			_mp.shownColumns.append(shown)
			if shown:
				noneShown = 0
			_mp._confDialog.newColumn(name, shown)

	if not noneShown:
		_mp._buildTable()

def readableName(model):
	if model.name:
		for char in model.name:
			if ord(char) < 32:
				break
		else:
			return model.name
	if isinstance(model, chimera.Molecule):
		return "unknown Molecule"
	if isinstance(model, chimera.MSMSModel):
		return "unknown MSMS surface"
	if isinstance(model, chimera.VRMLModel):
		return "unknown VRML object"
	return "unknown"

def inputPath(model):
	if not hasattr(model, 'openedAs') or '\n' in model.openedAs[0]:
		return readableName(model)
	path = model.openedAs[0]
	curdir = os.getcwd() + os.sep
	if path.startswith(curdir):
		return path[len(curdir):]
	return path

def getPhysicalChains(model):
	# return chains of physically connected residues as list of lists;
	# single-residue "chains" collated into first list
	from operator import add
	physical = [[]]
	seen = {}

	for root in model.roots(1):
		resAtoms = root.atom.residue.atoms
		numRootAtoms = root.size.numAtoms

		if numRootAtoms < len(resAtoms):
			# disconnected residue; continue only if this is
			# the largest fragment of the residue
			largestFrag = 1
			for atom in resAtoms:
				if atom.rootAtom == root.atom:
					continue
				if atom.molecule.rootForAtom(atom, 1).size\
						.numAtoms > numRootAtoms:
					largestFrag = 0
					break
			if not largestFrag:
				continue

		if numRootAtoms <= len(resAtoms):
			curPhysical = physical[0]
		else:
			curPhysical = []
			physical.append(curPhysical)
		for atom in model.traverseAtoms(root):
			res = atom.residue
			if seen.has_key(res):
				continue
			seen[res] = 1

			curPhysical.append(res)
	
	return physical

def nameColumn(m):
	if "trace" in readableName(m):
		if _mp and _mp._confDialog.showColorVar.get():
			bcolor = isinstance(m, chimera.Molecule) and m.color or None
			return readableName(m), bcolor
		return readableName(m)
	return ""
	
def _oslIdent(item):
	if isinstance(item, Group):
		models = item.models
		osls = [m.oslIdent() for m in models]
		from chimera.misc import oslModelCmp
		osls.sort(oslModelCmp)
		return u"%s\N{HORIZONTAL ELLIPSIS}%s" % (osls[0][1:], osls[-1][1:])
	return item.oslIdent()[1:]

from _surface import SurfaceModel
addColumns([
	('Name', 'text', nameColumn)
])
addColumns([
	('Note', 'text', lambda m: (hasattr(m, 'note') and m.note or '')),
	('Input file', 'text', inputPath)
], defaultShown=False)

class ModelPanel(ModelessDialog):
	title="SSETracer"
	buttons=('Run Tracer', 'Run Twister','Axis Comparison','Close')
	name="SSETracer"
	#help="UsersGuide/modelpanel.html"

	itemTableHelp = "click to select models;"\
			"\nright-hand action buttons work on selected models;"\
			"\ndouble-click to perform default action on model"\
			"\n(see 'Configure...' for default action info)"
	
	
			
	def __init__(self):
		self.outputQueue = Queue.Queue()
		self.runEnabled = True
		self.threadClass = None

		openModels = chimera.openModels.list()

		self.pdbCurrentPath = ""
		self.mrcCurrentPath = ""
		self.skeletonCurrentPath = ""
		self.stickCurrentPath = ""
		self.outputCurrentPath = ""
		self.pdbCurrentPathSSE = ""
		self.mrcCurrentPathSSE = ""
		self.skeletonCurrentPathSSE = ""

		for i in openModels:
			fileExt = os.path.splitext(i.openedAs[0])[1] #file extension
			if fileExt.lower() == ".pdb":
				self.pdbCurrentPath = i.openedAs[0]
			elif fileExt.lower() == ".mrc":
				if "skeleton" not in i.openedAs[0].lower():
					self.mrcCurrentPath = i.openedAs[0]
				else:
					self.skeletonCurrentPath = i.openedAs[0]

		#get path argument default from mrc path
		self.pdbPathSSE = Tkinter.StringVar()
		self.mrcPathSSE = Tkinter.StringVar()
		self.skeletonPathSSE = Tkinter.StringVar()
		self.pdbPathSSE.set(self.pdbCurrentPathSSE)
		self.mrcPathSSE.set(self.mrcCurrentPathSSE)
		self.skeletonPathSSE.set(self.skeletonCurrentPathSSE)
		
		self.pdbPath = Tkinter.StringVar()
		self.mrcPath = Tkinter.StringVar()
		self.skeletonPath = Tkinter.StringVar()
		self.stickPath = Tkinter.StringVar()
		self.outputPath = Tkinter.StringVar()
		self.pdbPath.set(self.pdbCurrentPath)
		self.mrcPath.set(self.mrcCurrentPath)
		self.skeletonPath.set(self.skeletonCurrentPath)
		self.stickPath.set(self.stickCurrentPath)
		self.outputPath.set(self.outputCurrentPath)
		#print "Current path: " + self.pdbCurrentPath

		self.pdb2Path = Tkinter.StringVar()
		self.mrc2Path = Tkinter.StringVar()
		self.pdb2Path.set(self.pdbCurrentPath)

		ModelessDialog.__init__(self)
		
	def fillInUI(self, parent):
		global _mp
		_mp = self
		self.parent = parent

		# model table
		self._getConfig()

		#not model panel
		from Tkinter import Entry, Label, Text, Frame, Button, Listbox, Checkbutton, Grid, Scrollbar
		
		self.label7SSE = Label(parent, text="SSETracer", font="-weight bold -underline 1")
		self.label7SSE.grid(row=0, column = 8, sticky='WS', pady=(12,2))
		self.label2SSE = Label(parent, text="MRC File:")
		self.label2SSE.grid(row=1, column = 7, sticky='E')
		self.mrcFileEntrySSE = Entry(parent, textvariable=self.mrcPathSSE)
		self.mrcFileEntrySSE.grid(row=1, column=8,sticky='EW')
		self.mrcFileButtonSSE = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.mrcPathSSE, [("MRC","*mrc")], ( (self.pdbPathSSE, ".pdb"), (self.skeletonPathSSE, "_skeleton.mrc")) ))
		self.mrcFileButtonSSE.grid(row=1, column=10,sticky='W')
		self.label3SSE = Label(parent, text="Skeleton File:")
		self.label3SSE.grid(row=2, column = 7, sticky='E')
		self.skeletonFileEntrySSE = Entry(parent, textvariable=self.skeletonPathSSE)
		self.skeletonFileEntrySSE.grid(row=2, column=8,sticky='EW')
		self.skeletonFileButtonSSE = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.skeletonPathSSE, [("MRC","*mrc")]))
		self.skeletonFileButtonSSE.grid(row=2, column=10,sticky='W')
		self.label4SSE = Label(parent, text="Threshold:")
		self.label4SSE.grid(row=3, column=7, sticky='E')
		self.thresholdEntrySSE= Entry(parent)
		self.thresholdEntrySSE.grid(row=3, column=8,sticky='EW')
		self.analysisSSE = Tkinter.IntVar()
		self.checkBox2SSE = Checkbutton(parent, variable=self.analysisSSE).grid(row=4, column =9, sticky='WS', pady=(13,0))
		self.label5SSE = Label(parent, text="Sensitivity analysis (using PDB below)")
		self.label5SSE.grid(row=4, column=8, sticky='WS', padx=(25,0), pady=(0,2))
		self.label1SSE = Label(parent, text="PDB File:")
		self.label1SSE.grid(row=5,column=7, sticky='E')
		self.pdbFileEntrySSE = Entry(parent, textvariable=self.pdbPathSSE)
		self.pdbFileEntrySSE.grid(row=5, column=8,sticky='EW')
		self.pdbFileButtonSSE = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.pdbPathSSE, [("PDB","*pdb")], ( (self.mrcPathSSE, ".mrc"), (self.skeletonPathSSE, "_skeleton.mrc") )  ))
		self.pdbFileButtonSSE.grid(row=5, column=10,sticky='W')
		
		#Strand Twister section
		self.label8SSE = Label(parent, text="Strand Twister", font="-weight bold -underline 1")
		self.label8SSE.grid(row=6, column=8, sticky='WS', pady=(28, 2))
		self.label1SSE = Label(parent, text="MRC File:")
		self.label1SSE.grid(row=7, column=7, sticky='E')
		self.mrc2FileEntrySSE = Entry(parent, textvariable=self.mrc2Path)
		self.mrc2FileEntrySSE.grid(row=7, column=8,sticky='EW')
		self.mrc2FileButtonSSE = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.mrc2Path, [("MRC","*mrc")], ((self.pdb2Path, ".pdb"),) ))
		self.mrc2FileButtonSSE.grid(row=7, column=10,sticky='W')
		self.label6SSE = Label(parent, text="Threshold:")
		self.label6SSE.grid(row=8, column=7, sticky='E')
		self.threshold2EntrySSE= Entry(parent)
		self.threshold2EntrySSE.grid(row=8, column=8,sticky='EW')
		self.threshold2EntrySSE.insert(0, '0') #default threshold
		self.label10SSE = Label(parent, text="PDB File:")
		self.label10SSE.grid(row=9, column=7, sticky='E')
		self.pdb2FileEntrySSE = Entry(parent, textvariable=self.pdb2Path)
		self.pdb2FileEntrySSE.grid(row=9, column=8,sticky='EW')
		self.pdb2FileButtonSSE = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.pdb2Path, [("PDB","*pdb")], ((self.mrc2Path, ".mrc"),) ))
		self.pdb2FileButtonSSE.grid(row=9, column=10,sticky='W')
		self.openCut = Tkinter.IntVar()
		
		self.label7 = Label(parent, text="Axis Comparison", font="-weight bold -underline 1")
		self.label7.grid(row=10, column = 8, sticky='WS', pady=(12,2))
		self.label12 = Label(parent, text="   ")
		self.label12.grid(row=11,column = 0, sticky='EW')
		self.label13 = Label(parent, text="   ")
		self.label13.grid(row=12,column = 0, sticky='EW')
		self.label14 = Label(parent, text="   ")
		self.label14.grid(row=13,column = 0, sticky='EW')
		self.label2 = Label(parent, text="PDB File:")
		self.label2.grid(row=11, column = 7, sticky='E')
		self.mrcFileEntry = Entry(parent, textvariable=self.mrcPath)
		self.mrcFileEntry.grid(row=11, column=8,sticky='EW')
		self.mrcFileButton = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.mrcPath, [("PDB","*pdb")] ))
		self.mrcFileButton.grid(row=11, column=10,sticky='W')
		self.label3 = Label(parent, text="First Helix(optional):")
		self.label3.grid(row=12, column = 7, sticky='E')
		self.skeletonFileEntry = Entry(parent, textvariable=self.skeletonPath)
		self.skeletonFileEntry.grid(row=12, column=8,sticky='EW')
		self.skeletonFileButton = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.skeletonPath, [("ALL","*")]))
		self.skeletonFileButton.grid(row=12, column=10,sticky='W')
		self.label4 = Label(parent, text="First Strand(optional):")
		self.label4.grid(row=13, column = 7, sticky='E')
		self.thresholdEntry = Entry(parent, textvariable=self.stickPath)
		self.thresholdEntry.grid(row=13, column=8,sticky='EW')
		self.thresholdButton = Button(parent,text="Browse", command=lambda: self.fileBrowse(self.stickPath, [("PDB","*pdb")]))
		self.thresholdButton.grid(row=13, column=10,sticky='W')
		self.label99 = Label(parent, text="Output(optional):")
		self.label99.grid(row=14, column = 7, sticky='E')
		self.outFileEntry = Entry(parent, textvariable=self.outputPath)
		self.outFileEntry.grid(row=14, column=8,sticky='EW')
		self.outFileButton = Button(parent,text="Browse", command=lambda: self.fileBrowse2(self.outputPath, [("TXT","*txt")] ))
		self.outFileButton.grid(row=14, column=10,sticky='W')
		self.label120 = Label(parent, text="              ")
		self.label120.grid(row=15,column = 0, sticky='EW')
		self.outputBox = Text(parent,height=16,width = 55)
		self.outputBox.grid(columnspan=9, row=16, sticky='E', pady=(10,0))#.rowconfigure(weight=1)
		self.scrollb = Scrollbar(parent, command=self.outputBox.yview)
		self.scrollb.grid(row=16, column=9, sticky='nsw', pady=(10,0))
		self.outputBox['yscrollcommand'] = self.scrollb.set
		
		# action buttons
		atf = self.allTitleFrame = Tkinter.Frame(self.parent)
		atf.grid(row=16, column=20, sticky='w')
		atf.grid_remove()
		self.commandLabel = Tkinter.Label(atf, text="Command")
		self.commandLabel.grid(row=0, column=1)
		Tkinter.Label(atf, text="Fav").grid(row=0, column=3, sticky='w')

		self.buttonScroll = Pmw.ScrolledFrame(self.parent,
							hscrollmode='none')
		self.buttonScroll.grid(row=16, column=20, sticky='nsew')
		self.favActionButtons = FavButtonBox(
			self.buttonScroll.interior(), orient='vertical', pady=0)
		self._shownActions = None
		self.allActionButtons = AllButtonBox(
			self.buttonScroll.interior(), self)
		self._favToggle = Pmw.RadioSelect(self.parent,
				command=self._favToggleCB,
				orient="horizontal", buttontype="radiobutton")
		self._favToggle.add("favorites")
		self._favToggle.add("all")
		self.favButtonsCreated = []
		self.allButtonsCreated = []

		self._addColumns()
		# add buttons from other extensions...
		self._addButtons()

		
		# add standard buttons
		def runClipping(models):
			import ModelClip.gui
			cd = chimera.dialogs.display(ModelClip.gui.ClipDialog.name)
			cd.setModel(models[0])
		addButton("close", openModels.close, moleculesOnly=False,
			balloon="close models")
		addButton("hide", lambda m, f='display', v=0,
			smf=setModelField: smf(m, f, v),
			moleculesOnly=False,
			balloon="hide selected models; undo with 'show'")
		def showRainbowDialog(models):
			from chimera import dialogs
			from rainbow import RainbowDialog
			if len(models) > 1:
				target = "models"
			else:
				target = "residues"
			dialogs.display(RainbowDialog.name).configure(
						models=models, target=target)
		addButton("select", selectCmd, moleculesOnly=False,
			balloon="incorporate models into graphics window"
			"\nselection using current selection mode"
			"\n(see graphics window Selection menu)")
		addButton("View All", viewerCmd, moleculesOnly=False,
			balloon="View all Models")
		from chainPicker import ChainPicker
		addButton("show", lambda m, f='display', v=1,
			smf=setModelField: smf(m, f, v),
			moleculesOnly=False, balloon="unhide selected models")
		def showTileDialog(models):
			from chimera import dialogs
			from EnsembleMatch.choose import TileStructuresCB
			dialogs.display(TileStructuresCB.name).configure(
								models=models)
		from transformDialog import TransformDialog
		from writePDBdialog import WritePDBdialog
		#from ksdsspDialog import KsdsspDialog
		
		self._favToggle.invoke("favorites")

		# add these last, since if they somehow fire before the
		# constructor is complete then an exception will occur
		chimera.triggers.addHandler('Model', self._fillTable, None)
		chimera.triggers.addHandler('OpenState', self._fillTable, None)
	
	def RunTracer(self):

		if self.runEnabled:

			self.outputBox.delete('1.0', Tkinter.END)

			threshold = self.thresholdEntrySSE.get()

			if self.mrcPathSSE and self.skeletonPathSSE and threshold:
				mrcPath = os.path.splitext(self.mrcPathSSE.get())[0] #file path and name without extension
				pdbPath = os.path.splitext(self.pdbPathSSE.get())[0] #file path and name without extension
				skeletonPath = os.path.splitext(self.skeletonPathSSE.get())[0] #file path and name without extension
				path = os.path.dirname(mrcPath)
				
				# Only use drive letter + colon if path points to drive root on Windows
				if platform.system() == 'Windows' and len(path) == 3:
					path = path[:2]

				if platform.system() == 'Linux':
					tracerPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + "tracer_v3_command"
				else:
					tracerPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + "tracer_v3_command.exe"

				tracerPath = '"' + tracerPath + '"'

				self.threadClass = TracerThread(tracerPath, path, mrcPath, skeletonPath, threshold, self.analysisSSE.get(), pdbPath, self.openCut.get(), self.outputQueue)
				self.threadClass.start()
				#buttons = ('Close')
				self.runEnabled = False
				#disable the actual button

				#waits to make sure thread.isAlive() returns true
				while True:
					try:
						if self.threadClass.isAlive():
							break
					except AttributeError:
						pass


				self.outputBox.after(1000, self.updateOutput)

			else:
				#change to dialog box
				self.outputBox.insert(Tkinter.END, "Incorrect input")

	def RunTwister(self):

		if self.runEnabled:

			self.outputBox.delete('1.0', Tkinter.END)

			threshold = self.threshold2EntrySSE.get()

			if self.mrc2Path and self.pdb2Path and threshold: #and mrcFile and newX and newY and newZ:
				mrcPath = os.path.splitext(self.mrc2Path.get())[0] #file path and name without extension
				pdbPath = os.path.splitext(self.pdb2Path.get())[0] #file path and name without extension
				#skeletonPath = os.path.splitext(self.skeletonPath.get())[0] #file path and name without extension
				path = os.path.dirname(mrcPath)

				if platform.system() == 'Linux':
					twisterPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + "strandtwister_v2_command"
				else:
					twisterPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + "strandtwister_v2_command.exe"

				twisterPath = '"' + twisterPath + '"'

				self.threadClass = TwisterThread(twisterPath, path, mrcPath, "skeletonPath", threshold, self.analysisSSE.get(), pdbPath, self.openCut.get(), self.outputQueue)
				self.threadClass.start()
				#buttons = ('Close')
				self.runEnabled = False
				#disable the actual button

				#waits to make sure thread.isAlive() returns true
				while True:
					try:
						if self.threadClass.isAlive():
							break
					except AttributeError:
						pass


				self.outputBox.after(1000, self.updateOutput)

			else:
				#change to dialog box
				self.outputBox.insert(Tkinter.END, "Incorrect input")

	#copy path from current entry to empty entry boxes in copyBoxes using the given extension from each tuple
	#copyBoxes contains one or more tuples in the format (DESTINATION, EXTENSION) where DESTINATION is a Var
	def copyToEmpty(self, currentBox, copyBoxes):
		for copyBox in copyBoxes:
				if not copyBox[0].get():
					filename = os.path.splitext(currentBox.get())[0] #file path and name without extension
					copyBox[0].set(filename + copyBox[1])

	#opens a file browser for the given types and then calls copyToEmpty to copy the filename with different extensions to the given locations
	def fileBrowse(self, currentBox, browseTypes, copyBoxes = []):
		from tkFileDialog import askopenfilename

		window = Tk()
		window.withdraw() # we don't want a full GUI, so keep the root window from appearing
		window.lift()
		window.attributes("-topmost", True)
		fullpath = os.path.normpath(
			askopenfilename(parent = window, filetypes = browseTypes))
		
		if fullpath and fullpath != ".":
			currentBox.set(fullpath)
			self.copyToEmpty(currentBox, copyBoxes)
			
	def fileBrowse2(self, currentBox, browseTypes, copyBoxes = []):
		from tkFileDialog import asksaveasfilename

		window = Tk()
		window.withdraw() # we don't want a full GUI, so keep the root window from appearing
		window.lift()
		window.attributes("-topmost", True)
		fullpath = asksaveasfilename(parent = window, filetypes = browseTypes)
		if fullpath and fullpath != ".":
			currentBox.set(fullpath)
			self.copyToEmpty(currentBox, copyBoxes)

	def AxisComparison(self):

		if self.runEnabled:

			self.outputBox.delete('1.0', Tkinter.END)

			#threshold = self.thresholdEntry.get()

			if self.mrcPath and self.skeletonPath:
				mrcPath = os.path.splitext(self.mrcPath.get())[0] #file path and name without extension
				pdbPath = os.path.splitext(self.pdbPath.get())[0] #file path and name without extension
				skeletonPath = os.path.splitext(self.skeletonPath.get())[0] #file path and name without extension
				stickPath = os.path.splitext(self.stickPath.get())[0]
				outPath = os.path.splitext(self.outputPath.get())[0]
				path = os.path.dirname(mrcPath)

				if platform.system() == 'Linux':
					tracerPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + "leastsquare"
				else:
					tracerPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + "leastsquare.exe"
					
				tracerPath = '"' + tracerPath + '"'
				
				proteinNameOffset = 1	
				stringHolder = mrcPath[-proteinNameOffset:]
				while ("\\" not in stringHolder):
					proteinNameOffset+=1
					stringHolder = mrcPath[-proteinNameOffset:]
				proteinNameOffset -= 1	
				
				if(outPath[-4:] != ".txt") and (outPath != "") :
					outPath = outPath + ".txt"
				# close all currently open models
				for e in chimera.openModels.list(all=True):
					chimera.openModels.close(e)
				#delete old files
				if os.path.exists(mrcPath[:-proteinNameOffset] + "output/"):
					shutil.rmtree(mrcPath[:-proteinNameOffset] + "output/") 
				if not os.path.exists(mrcPath[:-proteinNameOffset] + "/output/"):
					os.makedirs(mrcPath[:-proteinNameOffset] + "/output/")
					
				self.threadClass = DistanceCompareThread(tracerPath, path, mrcPath, skeletonPath, stickPath, pdbPath, self.outputQueue, outPath)
				
				self.threadClass.start()
				#buttons = ('Close')
				self.runEnabled = False
				#disable the actual button

				#waits to make sure thread.isAlive() returns true
				while True:
					try:
						if self.threadClass.isAlive():
							break
					except AttributeError:
						pass


				self.outputBox.after(1000, self.updateOutput)

			else:
				#change to dialog box
				self.outputBox.insert(Tkinter.END, "Incorrect input")
		
		#wait for c++ program to create files
		while self.threadClass.isAlive():
			True
		#open files
		it = 0
		x = -1
		if skeletonPath != "" :
			while(skeletonPath[len(skeletonPath)-it-1] != "1" and skeletonPath[len(skeletonPath)-it-1] != "0"):
				if(len(skeletonPath)-it-2 > 0):
					it = it+1
				else:
					break;
			if(skeletonPath[len(skeletonPath)-it-1] == "1"):
				x = 1
				if(os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/traceHelix0.pdb")):
					os.remove(mrcPath[:-proteinNameOffset] +  "/output/traceHelix0.pdb")
			if(skeletonPath[len(skeletonPath)-it-1] == "0"):
				x = 0
		if stickPath != "" :
			while(stickPath[len(stickPath)-it-1] != "1" and stickPath[len(stickPath)-it-1] != "0"):
				it = it+1
			if(stickPath[len(stickPath)-it-1] == "1"):
				x = 1
				#if(os.path.isfile(mrcPath[:-4] +  "/output/traceHelix0.pdb")):
				#	os.remove(mrcPath[:-4] +  "/output/traceHelix0.pdb")
			if(stickPath[len(stickPath)-it-1] == "0"):
				x = 0
		
		x = 1
		if (os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/traceHelix0" + ".pdb")):
			chimera.openModels.open(mrcPath[:-proteinNameOffset] + "/output/traceHelix0" + ".pdb")
		#time.sleep(.5)
		count = x
		if (os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/traceHelix" + str(count) + ".pdb")):
			while(True):
				while(os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/traceHelix" + str(count) + ".pdb")):
					chimera.openModels.open(mrcPath[:-proteinNameOffset] + "/output/traceHelix" + str(count) + ".pdb")
					count = count + 1	
				if(count != x):
					break
				else:
					time.sleep(.25)
			
		variable = 0
		count = x
		
		if (os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/traceSheet0" + ".pdb")):
			chimera.openModels.open(mrcPath[:-proteinNameOffset] + "/output/traceSheet0" + ".pdb")
		while(True):
			while(os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/traceSheet" + str(count) + ".pdb")):
				chimera.openModels.open(mrcPath[:-proteinNameOffset] + "/output/traceSheet" + str(count) + ".pdb")
				count = count + 1	
				variable = 0
			if(time > 4):
				break
			if(count != x):
				break
			else:
				variable = variable + 1
				time.sleep(.25)		
				
	
		count = 1
		while(os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/trueHelix" + str(count) + ".pdb")):
			chimera.openModels.open(mrcPath[:-proteinNameOffset] + "/output/trueHelix" + str(count) + ".pdb")
			count = count + 1
		while(os.path.isfile(mrcPath[:-proteinNameOffset] +  "/output/trueSheet" + str(count) + ".pdb")):
			chimera.openModels.open(mrcPath[:-proteinNameOffset] + "/output/trueSheet" + str(count) + ".pdb")
			count = count + 1
		
		if (os.path.isfile(skeletonPath + ".mrc")):
			chimera.openModels.open(skeletonPath + ".mrc")
			for v in openModels.list(modelTypes=[Volume]):
				v.initialize_thresholds ( self,  True )
				v.set_parameters(surface_colors = [(0.7, 0.7, 0.7, 0.5)], surface_levels = [.1])
				v.show()
			
		chimera.openModels.open(mrcPath + ".pdb")
		from chimera import viewer
		viewer.viewAll()
		
		import Midas
		Midas.color('#d2d2b4b48c8c', '@/element=C')

	def updateOutput(self):

		try:
			#self.outputBox.insert(Tkinter.END, "try block")
			self.outputBox.insert(Tkinter.END, self.outputQueue.get(block=False))
			self.outputBox.yview(Tkinter.MOVETO, 1.0)
		except Queue.Empty:
			pass

		#self.outputBox.insert(Tkinter.END, "test")
		try:
			if not self.runEnabled:
				self.outputBox.after(50, self.updateOutput)

			if not self.threadClass.isAlive(): #make sure 
			 	self.outputBox.after(3000, self.enableButton)
		except AttributeError:
			pass

	def enableButton(self):
		while True:
			try:
				self.outputBox.insert(Tkinter.END, self.outputQueue.get(block=False))
				self.outputBox.yview(Tkinter.MOVETO, 1.0)
			except Queue.Empty:
				break

		self.runEnabled = True
		#enable the actual button

		#self.update_idletasks()

		#call updateOutput until process finished
		#disable run button until finished
		#add stop button
		
	def Configure(self):
		"""configure action buttons"""
		self._confDialog.enter()

	def see(self, buttonName):
		pass # I don't think anything calls this

	def selected(self, moleculesOnly=False, groupsOkay=False):
		"""Return a list of the selected models"""

		selected = []
		for ii in self.itemTable.hlist.info_selection():
			item = self.items[int(ii)]
			if groupsOkay:
				selected.append(item)
				continue
			models = _getModels(item)
			if moleculesOnly:
				models = [m for m in models if isinstance(m, Molecule)]
			selected.extend(models)
		return selected

	def selectionChange(self, models, extend=False, priorSelection=None):
		"""set (or extend) the selection to contain the given models
		
		   'models' can be Models or oslIdents"""

		# may have to ungroup groups if they are partially selected
		newSelected = []
		breakGroups = []
		if models:
			testSet = set()
			for m in models:
				testSet.update(_getModels(m))
			if isinstance(models[0], basestring):
				# OSL ident
				for i, item in enumerate(self.items):
					val = item.oslIdent()
					if isinstance(val, set):
						if testSet & val:
							breakGroups.append(item)
					elif val in testSet:
						newSelected.append(i)
			else:
				for i, item in enumerate(self.items):
					val = set(_getModels(item))
					if val & testSet:
						if val & testSet < val:
							breakGroups.append(item)
						else:
							newSelected.append(i)
		if breakGroups:
			if priorSelection is None:
				priorSelection = self.selected(groupsOkay=True)
			for group in breakGroups:
				self.items.remove(group)
				self.items.extend(group.components)
				if group in priorSelection:
					priorSelection.remove(group)
					priorSelection.extend(group.components)
			self.selectionChange(models, extend=extend,
					priorSelection=priorSelection)
			return
		if priorSelection is not None:
			self._fillTable(fromScratch=True, selected=priorSelection)
		if not extend:
			self.itemTable.hlist.selection_clear()

		for i in newSelected:
			self.itemTable.hlist.selection_set(i)
		self._selChangeCB()

	def _addButtons(self):
		"""Add buttons to interface that were requested before
		   panel was created.
		"""

		for name, info in _buttonInfo.items():
			balloon, defaultFavorite = info[-3:-1]
			self._confDialog.newButton(name, balloon=balloon,
						defaultFavorite=defaultFavorite)
			self._showButton(name)

	def _addColumns(self):
		"""Process column information"""
		self.shownColumns = []

		for i in range(len(_columnNames)):
			name = _columnNames[i]
			if name == "Note":
				shown = False
				for m in openModels.list():
					if hasattr(m, 'note') and m.note:
						shown = True
						break
			else:
				try:
					shown = self._confDialog.prefs[
							'shownColumns'][name]
				except KeyError:
					shown = _defaultShows[i]
			self.shownColumns.append(shown)
			self._confDialog.newColumn(name, shown)
		self._buildTable()

	def _buttonParams(self, name):
		callback, minModels, maxModels, moleculesOnly, balloon, \
					defaultFavorite, groupsOkay = _buttonInfo[name]
		kw = {}
		state = 'normal'
		if self._shouldDisable(minModels, maxModels, moleculesOnly):
			state = 'disabled'
		kw['state'] = state
		kw['pady'] = 0
		# if you click a button fast enough, you can get around it's
		# upcoming disabling...
		def cmd(cb=callback, s=self, mo=moleculesOnly, minm=minModels,
				maxm=maxModels, go=groupsOkay):
			if not s._shouldDisable(minm, maxm, mo):
				cb(s.selected(moleculesOnly=mo, groupsOkay=go))
		kw['command'] = cmd
		return kw, balloon, defaultFavorite

	def _buildTable(self):
		if hasattr(self, 'itemTable'):
			# can't dynamically add columns to Tix widget;
			# destroy and recreate
			selected = self.selected()
			self.itemTable.grid_forget()
			self.itemTable.destroy()
		else:
			selected = None

		w, h = self._confDialog.prefs['table w/h']
		inch = self.parent.winfo_fpixels("1i")
		self.itemTable = Tix.ScrolledHList(self.parent,
			width=150,
			height=50,
			options="""hlist.columns %d
			hlist.header 1
			hlist.selectMode extended
			hlist.indicator 0"""
			% len(filter(lambda s: s == 1, self.shownColumns)))
		help.register(self.itemTable, balloon=self.itemTableHelp)
		self.itemTable.hlist.config(browsecmd=self._selChange,
							command=self._dblClick)
		self.textStyle = Tix.DisplayStyle("text",
				refwindow=self.itemTable)
		# get a style for checkbutton columns...
		self.checkButtonStyle = Tix.DisplayStyle("window",
				refwindow=self.itemTable, anchor="center")
		self.colorWellStyle = Tix.DisplayStyle("window",
				refwindow=self.itemTable, anchor="center")
		colNum = 0
		self.columnMap = []
		showFullTitles = False
		last = self._confDialog.prefs["lastUse"]
		from time import time
		now = self._confDialog.prefs["lastUse"] = time()
		if last is None or now - last > 777700: # about 3 months
			showFullTitles = True
		for index in range(len(_columnNames)):
			if not self.shownColumns[index]:
				continue
			self.columnMap.append(index)
			text = _columnNames[index]
			if _valueTypes[index] == 'toggle' \
			and not showFullTitles:
				text = text[:1]
			self.itemTable.hlist.header_create(colNum,
						itemtype='text', text=text)
			colNum = colNum + 1
			
		self.parent.columnconfigure(20, weight=1)
		self.parent.rowconfigure(20, weight=1)
		self.itemTable.grid(row=16, column=10, sticky='nsew',
								rowspan=1)
		self._fillTable(selected=selected, fromScratch=1)
		self.itemTable.bind("<Configure>", self._rememberSize, add=True)

	def _dblClick(self, item):
		"""user has double-clicked on model table entry"""

		# if the state of the action buttons is due to change,
		# execute that change before calling the double-click routine
		if hasattr(self, '_selChangeIdle') and self._selChangeIdle:
			self.parent.after_cancel(self._selChangeIdle)
			self._selChangeCB()

		self._confDialog.dblClick()

	def _favButton(self, name, fav):
		names = self.favButtonsCreated
		actionButtons = self.favActionButtons
		if fav:
			names.append(name)
			names.sort(lambda a, b: cmp(a.lower(), b.lower()))
			kw, balloon, defaultFavorite = self._buttonParams(name)
			index = names.index(name)
			if index == len(names)-1:
				addFunc = actionButtons.add
			else:
				addFunc = actionButtons.insert
				kw['beforeComponent'] = names[index+1]

			but = addFunc(name, **kw)
			but.config(default='disabled')
			if balloon:
				help.register(but, balloon=balloon)
		else:
			names.remove(name)
			actionButtons.delete(name)

	def _favToggleCB(self, label):
		if label == "favorites":
			self._showActions(self.favActionButtons)
		else:
			self._showActions(self.allActionButtons)

	def _fillTable(self, *triggerArgs, **kw):
		if len(triggerArgs) > 0:
			if triggerArgs[0] == 'OpenState':
				if 'active change' not in triggerArgs[-1].reasons:
					return
			elif triggerArgs[0] == 'Model':
				global _groupNameCache
				_groupNameCache.clear()
		hlist = self.itemTable.hlist
		defaultable = False
		if kw.get('selected', None) != None:
			selected = kw['selected']
		else:
			selected = self.selected(groupsOkay=True)
			defaultable = True
		rebuild = False
		curModels = set(openModels.list())
		if not hasattr(self, 'items'):
			global _groups
			self.items = _groups[:]
			_groups[:] = []
			rebuild = True
		elif kw.get('fromScratch', False):
			rebuild = True
		else:
			prevModels = set()
			for item in self.items:
				prevModels.update(_getModels(item))
			if curModels != prevModels:
				rebuild = True

		if rebuild:
			newItems = []
			for item in self.items:
				if not isinstance(item, Group):
					continue
				item.update()
				if len(item.models) > 1:
					newItems.append(item)
					curModels.difference_update(item.models)
			self.items = newItems + list(curModels)
			self.items.sort(self._itemSort)
			self._prevValues = {}
			hlist.delete_all()
			vf = _valueFuncs[self.columnMap[0]]
			for i, item in enumerate(self.items):
				hlist.add(i, **self._hlistKw(item, 0))
				self._prevValues[(i, 0)] = vf(item)
			for ci in range(1, len(self.columnMap)):
				vf = _valueFuncs[self.columnMap[ci]]
				for i, item in enumerate(self.items):
					hlist.item_create(i, ci, **self._hlistKw(item, ci))
					self._prevValues[(i, ci)] = vf(item)
		else:
			for ci in range(len(self.columnMap)):
				vf = _valueFuncs[self.columnMap[ci]]
				for i, item in enumerate(self.items):
					curVal = vf(item)
					prevVal = self._prevValues[(i, ci)]
					if isinstance(curVal, tuple):
						for vi in range(len(curVal)):
							valItem = curVal[vi]
							if callable(valItem) and not isinstance(valItem,
									GroupAttr):
								continue
							pv = prevVal[vi]
							if type(valItem) != type(pv) or valItem != pv:
								break
						else:
							# equal
							continue
					elif curVal == prevVal:
						continue
					self._prevValues[(i, ci)] = curVal
					hlist.item_configure(i, ci, **self._hlistKw(item, ci))
		# if only one item, select it
		if defaultable and len(self.items) == 1:
			selected = self.items
		for item in selected:
			if item not in self.items:
				continue
			hlist.selection_set(self.items.index(item))
		self._selChange(None)

	def _getConfig(self):
		"""retrieve configuration preferences"""

		# set up configuration dialog
		from confDialog import ConfDialog
		self._confDialog = ConfDialog(self)
		self._confDialog.Close()

	def _hlistKw(self, item, colNum):
		vt = _valueTypes[self.columnMap[colNum]]
		vf = _valueFuncs[self.columnMap[colNum]]
		kw = {'itemtype': vt}
		txt = None
		img = None
		val = vf(item)
		if isinstance(val, set) and vt not in ['toggle', 'well']:
			return {}
		from Group import GroupAttr
		if vt == 'text':
			txt = val
			if isinstance(txt, GroupAttr):
				testable = list(txt.vals)[0]
			else:
				testable = txt
			if not isinstance(testable, basestring):
				txt, bcolor = txt
				if bcolor is not None:
					if not isinstance(bcolor, basestring):
						if hasattr(bcolor, 'rgba'):
							rgba = bcolor.rgba()
						else:
							rgba = bcolor
						from CGLtk.color import rgba2tk
						bcolor = rgba2tk(rgba)
						fcolor = CGLtk.textForeground(
							bcolor, self.itemTable)
					kw['style'] = Tix.DisplayStyle("text",
						refwindow=self.itemTable,
						background=bcolor,
						foreground=fcolor,
						selectforeground=bcolor)
			else:
				kw['style'] = self.textStyle
		elif vt == 'image':
			img = val
		elif vt == 'imagetext':
			img, txt = val
		elif vt == 'toggle':
			kw['itemtype'] = 'window'
			truth, cb = val
			togKw = {'command':
				# avoid holding references to model
				lambda cb=cb, i=self.items.index(item),
					nt=isinstance(truth, GroupAttr) or not truth:
					cb(self.items[i], nt),
				'indicatoron': 0,
				'borderwidth': 0}
			if isinstance(truth, GroupAttr):
				togKw['image'] = self.itemTable.tk.call(
					'tix', 'getimage', 'ck_onoff_37')
			elif truth:
				togKw['image'] = self.itemTable.tk.call(
					'tix', 'getimage', 'ck_on')
			else:
				togKw['image'] = self.itemTable.tk.call(
					'tix', 'getimage', 'ck_off')
			toggle = Tkinter.Checkbutton(
						self.itemTable.hlist, **togKw)
			kw['window'] = toggle
			kw['style'] = self.checkButtonStyle
		elif vt == 'well':
			color, noneOkay, alphaOkay, cb = val
			if color is False:
				kw['itemtype'] = 'text'
				txt = ""
			else:
				kw['itemtype'] = 'window'
				if isinstance(color, chimera.MaterialColor):
					color = color.rgba()
				from weakref import proxy
				def wellCB(clr, cb=cb, mdl=proxy(item)):
					if clr is not None:
						clr = chimera.MaterialColor(
									*clr)
					cb(mdl, clr)
				from CGLtk.color.ColorWell import ColorWell
				kw['window'] = ColorWell(self.itemTable.hlist,
					color, callback=wellCB,
					multiple=isinstance(color, GroupAttr),
					width=18, height=18,
					noneOkay=noneOkay, wantAlpha=alphaOkay)
				kw['style'] = self.colorWellStyle
		else:
			raise ValueError("Unknown column type: '%s'" % vt)
		
		if txt != None:
			kw['text'] = unicode(txt)
		if img != None:
			kw['image'] = self.itemTable.tk.call(
							'tix', 'getimage', img)
		return kw
	
	def _itemSort(self, i1, i2):
		def getVal(vals):
			if isinstance(vals, set):
				return min(vals)
			return vals
		id1 = getVal(i1.id)
		id2 = getVal(i2.id)
		if id1 < id2:
			return -1
		if id1 > id2:
			return 1
		subid1 = getVal(i1.subid)
		subid2 = getVal(i2.subid)
		if subid1 < subid2:
			return -1
		if subid1 > subid2:
			return 1
		return 0

	def _rememberSize(self, event):
		w, h = event.width, event.height
		if min(w,h) < 20:
			return
		inch = self.parent.winfo_fpixels("1i")
		self._confDialog.prefs["table w/h"] = (w/inch, h/inch)

	def _selChange(self, item):
		# slow browse callback interferes with double-click detection,
		# so delay callback enough to allow most double-clicks to work
		if hasattr(self, '_selChangeIdle') and self._selChangeIdle:
			self.parent.after_cancel(self._selChangeIdle)
		self._selChangeIdle = self.parent.after(300, self._selChangeCB)

	def _selChangeCB(self):
		numSel = len(self.itemTable.hlist.info_selection())
		allButtons = _buttonInfo.keys()
		favs = self._confDialog.prefs["favorites"]
		for buttons, actionButtons in [
				([b for b in allButtons if favs[b]], self.favActionButtons),
				(allButtons, self.allActionButtons)]:
			for b in buttons:
				state = 'normal'
				callback, minModels, maxModels, moleculesOnly, \
					balloon, defaultFavorite, groupsOkay \
							= _buttonInfo[b]
				if self._shouldDisable(minModels, maxModels,
								moleculesOnly):
					state = 'disabled'
				actionButtons.button(b).config(state=state)
		self._selChangeIdle = None

	def _shouldDisable(self, minModels, maxModels, moleculesOnly):
		if moleculesOnly:
			numSel = len(self.selected(moleculesOnly=True))
		else:
			numSel = len(self.itemTable.hlist.info_selection())
		if minModels != None and numSel < minModels \
		or maxModels != None and numSel > maxModels:
			return 1
		return 0

	def _showActions(self, actionButtons):
		if actionButtons == self._shownActions:
			return
		if self._shownActions:
			if actionButtons == self.favActionButtons:
				self.allActionButtons.grid_remove()
				self.allTitleFrame.grid_remove()
			else:
				self.favActionButtons.grid_remove()
				self.allTitleFrame.grid()
				bw = actionButtons.buttonWidth()
				lw = self.commandLabel.winfo_reqwidth()
				pad = (bw-lw) /2.0
				self.allTitleFrame.columnconfigure(0, minsize=pad)
				self.allTitleFrame.columnconfigure(2, minsize=pad)
		actionButtons.grid()
		self.buttonScroll.component('clipper').configure(
						width=actionButtons.clipWidth()+2, height='2i')
		if chimera.tkgui.windowSystem == 'aqua':
			# work around bug where Aqua would behave as if the
			# scroller was to the right of its actual position
			# once the clipper was narrowed
			def later(tl = self.allTitleFrame.winfo_toplevel()):
				tl.wm_geometry(tl.wm_geometry())
			self.allTitleFrame.after(100, later)
		self.buttonScroll.yview(mode='moveto', value=0.0)
		self._shownActions = actionButtons

	def _showButton(self, name):
		kw, balloon, defaultFavorite = self._buttonParams(name)
		favPrefs = self._confDialog.prefs['favorites']
		names = self.allButtonsCreated
		actionButtons = self.allActionButtons
		names.append(name)
		names.sort(lambda a, b: cmp(a.lower(), b.lower()))
		index = names.index(name)
		if index == len(names)-1:
			addFunc = actionButtons.add
		else:
			addFunc = actionButtons.insert
			kw['beforeComponent'] = names[index+1]

		but = addFunc(name, **kw)
		but.config(default='disabled')
		if balloon:
			help.register(but, balloon=balloon)
		if favPrefs.get(name, defaultFavorite):
			self._favButton(name, True)

if not chimera.nogui:
	class FavButtonBox(Pmw.ButtonBox):
		def __init__(self, *args, **kw):
			Pmw.ButtonBox.__init__(self, *args, **kw)

		def clipWidth(self):
			maxWidth = 0
			for i in range(self.numbuttons()):
				w = self.button(i).winfo_reqwidth()
				if w > maxWidth:
					maxWidth = w
			return maxWidth

		buttonWidth = clipWidth

	class AllButtonBox(Tkinter.Frame):
		def __init__(self, master, modelPanel):
			Tkinter.Frame.__init__(self, master)
			self.__rowInfo = {}
			self.__modelPanel = modelPanel
			self.__maxButtonWidth = 0
			

		def add(self, name, **buttonKw):
			return self.__addButton(name, len(self.__rowInfo), **buttonKw)

		def button(self, name):
			return self.__rowInfo[name][2]

		def clipWidth(self):
			if self.__rowInfo:
				return self.__maxButtonWidth + self.__checkWidth
			return 0

		def buttonWidth(self):
			return self.__maxButtonWidth

		def insert(self, name, beforeComponent=None, **buttonKw):
			at = self.__rowInfo[beforeComponent][0]
			for bname, info in self.__rowInfo.items():
				row, frame, button, check = info
				if row >= at:
					frame.grid_forget()
					frame.grid(row=row+1, column=0)
					self.__rowInfo[bname][0] += 1
			return self.__addButton(name, at, **buttonKw)

		def __addButton(self, name, row, **buttonKw):
			f = Tkinter.Frame(self)
			f.grid(row=row, column=0)
			b = Tkinter.Button(f, text=name, **buttonKw)
			b.grid(row=0, column=0, sticky="ew")
			bw = b.winfo_reqwidth()
			if bw > self.__maxButtonWidth:
				for info in self.__rowInfo.values():
					info[1].columnconfigure(0, minsize=bw)
				self.__maxButtonWidth = bw
			else:
				f.columnconfigure(0, minsize=self.__maxButtonWidth)
			chkKw = {'command': lambda nm=name, s=self: s.__changeFav(nm),
				'indicatoron': 0, 'pady': 0,
				'borderwidth': 0}
			isFav = self.__modelPanel._confDialog.prefs['favorites'][name]
			if isFav:
				chkKw['image'] = f.tk.call('tix', 'getimage', 'ck_on')
			else:
				chkKw['image'] = f.tk.call('tix', 'getimage', 'ck_off')
			if chimera.tkgui.windowSystem == 'aqua':
				chkKw['relief'] = 'flat'
				chkKw['bd'] = 0
			check = Tkinter.Checkbutton(f, **chkKw)
			self.__checkWidth = check.winfo_reqwidth()
			check.grid(row=0, column=1)
			self.__rowInfo[name] = [row, f, b, check]
			return b

		def __changeFav(self, name):
			favs = self.__modelPanel._confDialog.prefs['favorites']
			fav = favs[name]
			self.__modelPanel._favButton(name, not fav)
			favsCopy = favs.copy()
			favsCopy[name] = not fav
			self.__modelPanel._confDialog.prefs['favorites'] = favsCopy
			chk = self.__rowInfo[name][-1]
			if fav:
				chk.configure(image=chk.tk.call('tix', 'getimage', 'ck_off'))
			else:
				chk.configure(image=chk.tk.call('tix', 'getimage', 'ck_on'))

from chimera import dialogs
dialogs.register(ModelPanel.name, ModelPanel)

def _setAttr(m, field, value, openState=0):
	if openState:
		setattr(m.openState, field, value)
	else:
		setattr(m, field, value)

# functions used in model panel button; could be called directly also
def setModelField(models, field, value, openState=0):
	for m in models:
		_setAttr(m, field, value, openState)

def setModelFieldOnly(models, field, onVal=1, offVal=0, openState=0):
	# turn off first, then on, so that models not in the models list
	# that nonetheless have shared openStates get the 'on' value
	for m in openModels.list():
		_setAttr(m, field, offVal, openState)
	for m in models:
		_setAttr(m, field, onVal, openState)

def toggleModelField(models, field, onVal=1, offVal=0, openState=0):
	openStates = {}
	for m in models:
		if openState:
			openStates[m.openState] = 1
			continue
		if curval == onVal:
			_setAttr(m, field, offVal, openState)
		else:
			_setAttr(m, field, onVal, openState)
	for os in openStates.keys():
		if getattr(os, field) == onVal:
			setattr(os, field, offVal)
		else:
			setattr(os, field, onVal)

_prevActivities = None
def activateAllCmd():
	"""Activate all models.  Restore previous activities if called again."""

	if _mp:
		favPrefs = _mp._confDialog.prefs['favorites']
		if favPrefs['activate all']:
			actionButtons = _mp.favActionButtons
		else:
			actionButtons = _mp.allActionButtons
	global _prevActivities
	if _prevActivities:
		for m in openModels.list():
			if _prevActivities.has_key(m.openState):
				m.openState.active = _prevActivities[
								m.openState]
		_prevActivities = None
		if _mp:
			butText = 'activate all'
	else:
		_prevActivities = {}
		for m in openModels.list():
			if _prevActivities.has_key(m.openState):
				continue
			_prevActivities[m.openState] = m.openState.active
			m.openState.active = 1
		if _mp:
			butText = 'restore activities'
	if _mp:
		actionButtons.component('activate all').config(text=butText)

_attrInspectors = {}
_headers = {}
_seqInspectors = {}
_inspectors = [_attrInspectors, _headers] # _seqInspectors is per chain

def _checkTrigger():
	global _modelTrigger
	for inspDict in _inspectors:
		if len(inspDict) > 0:
			break
	else:
		# should be no trigger active; start one
		_modelTrigger = chimera.triggers.addHandler(
						'Model', _modelTriggerCB, None)

def attributesCmd(models):
	global _attrInspectors
	_checkTrigger()
	for model in models:
		if not _attrInspectors.has_key(model):
			from modelInspector import ModelInspector
			_attrInspectors[model] = ModelInspector(model)
		_attrInspectors[model].enter()

def seqCmd(items):
	global _seqInspectors
	todo = []
	for item in items[:]:
		if not _seqInspectors.has_key(item):
			from chimera.Sequence import StructureSequence
			if isinstance(item, StructureSequence):
				copySeq = StructureSequence.__copy__(item)
				copySeq.name = item.fullName()
				_addSeqInspector(item, mavSeq=copySeq)
			else:
				seqs = item.sequences()
				if seqs:
					todo.extend(seqs)
				else:
					items.remove(item)
				continue
		_seqInspectors[item].enter()
	if todo:
		if len(todo) > 1:
			from seqPanel import SeqPickerDialog
			from chimera import dialogs
			d = dialogs.display(SeqPickerDialog.name)
			d.molListBox.setvalue(todo)
		else:
			seqCmd(todo)
	else:
		# if handed only sequences...
		return [_seqInspectors[item] for item in items]

def _addSeqInspector(item, mavSeq=None, mav=None):
	global _seqInspectors, _saveSessionTrigger
	if not _seqInspectors:
		from SimpleSession import SAVE_SESSION
		_saveSessionTrigger = chimera.triggers.addHandler(
				SAVE_SESSION, _saveSessionCB, None)
	trigMav = []
	hid = item.triggers.addHandler(item.TRIG_DELETE,
			lambda tn, md, td: md[0].Quit(), trigMav)
	def quitCB(mav, i=item, h=hid):
		del _seqInspectors[i]
		i.triggers.deleteHandler(i.TRIG_DELETE, h)
		if not _seqInspectors:
			from SimpleSession import SAVE_SESSION
			chimera.triggers.deleteHandler(SAVE_SESSION, _saveSessionTrigger)
	if mav:
		mav.quitCB = quitCB
	else:
		from MultAlignViewer.MAViewer import MAViewer
		if mavSeq.descriptiveName:
			title = "chain %s: %s" % (mavSeq.chainID, mavSeq.descriptiveName)
		else:
			title = mavSeq.name
		mav = MAViewer([mavSeq], title=title, quitCB=quitCB,
						autoAssociate=None, sessionSave=False)
	_seqInspectors[item] = mav
	trigMav.append(mav)

def _modelTriggerCB(trigName, myArg, modelsChanges):
	for model in modelsChanges.deleted:
		for inspDict in _inspectors:
			if inspDict.has_key(model):
				_deleteInspector(model, inspDict)

def _saveSessionCB(trigName, myArg, session):
	from SimpleSession import sessionID, sesRepr
	info = []
	for seq, mav in _seqInspectors.items():
		info.append((seq.name, sessionID(seq.molecule),
			[seq.saveInfo() for seq in mav.seqs], mav.saveInfo()))
	print>>session, """
try:
	from ModelPanel import restoreSeqInspectors
	restoreSeqInspectors(%s)
except:
	reportRestoreError("Error restoring sequence viewers")
""" % sesRepr(info)

def restoreSeqInspectors(info):
	global _seqInspectors
	from SimpleSession import idLookup
	for seqName, molID, seqsInfo, mavInfo in info:
		mol = idLookup(molID)
		for seq in mol.sequences():
			if seq.name == seqName:
				break
		else:
			continue
		if seq in _seqInspectors:
			continue

		from MultAlignViewer.MAViewer import restoreMAV
		from chimera.Sequence import restoreSequence
		mav = restoreMAV([restoreSequence(seqInfo)
						for seqInfo in seqsInfo], mavInfo)
		_addSeqInspector(seq, mav=mav)

def _deleteInspector(model, dict):
	inspector = dict[model]
	del dict[model]
	for inspDict in _inspectors:
		if len(inspDict) > 0:
			break
	else:
		# no inspectors; drop triggers
		chimera.triggers.deleteHandler('Model', _modelTrigger)
	if hasattr(inspector, 'destroy'):
		inspector.destroy()
	else:
		inspector._toplevel.destroy()

def backboneCmd(models, resTrace=1):
	from chimera.misc import displayResPart
	for m in models:
		if not hasattr(m, 'residues'):
			continue
		if resTrace:
			displayResPart(m.residues, trace=1)
		else:
			displayResPart(m.residues, backbone=1)

def focusCmd(models):
	from chimera import openModels, viewer, update
	shown = {}
	for m in openModels.list():
		shown[m] = m.display
		if m in models:
			m.display = 1
		else:
			m.display = 0
	update.checkForChanges()
	viewer.viewAll()
	if chimera.openModels.cofrMethod != chimera.OpenModels.Independent:
		openModels.cofrMethod = openModels.CenterOfView
		viewer.clipping = True

	for m,disp in shown.items():
		m.display = disp
	update.checkForChanges()

_groups = []
_groupNameCache = {}
def groupCmd(items, name=None):
	from Group import Group
	removeGroups = []
	addItems = []
	newGroup = None
	if len(items) == 1:
		item = items[0]
		if isinstance(item, Group):
			removeGroups.append(item)
			addItems.extend(item.components)
			sel = item.components
		else:
			from chimera import UserError
			raise UserError("Cannot group a single model.")
	else:
		models = []
		for item in items:
			if isinstance(item, Group):
				removeGroups.append(item)
				models.extend(item.models)
			else:
				models.append(item)
		global _groupNameCache
		models.sort()
		key = tuple([id(m) for m in models])
		if name is None:
			if key in _groupNameCache:
				name = _groupNameCache[key]
			else:
				name = getGroupName(items)
		_groupNameCache[key] = name
		newGroup = Group(items, name)
		addItems.append(newGroup)
		sel = addItems
	if _mp:
		selected = _mp.selected()
		for rg in removeGroups:
			_mp.items.remove(rg)
		_mp.items.extend(addItems)
		_mp._fillTable(fromScratch=True, selected=selected)
		_mp.selectionChange(sel)
	else:
		if addItems:
			addGroups = [ai for ai in addItems if isinstance(ai, Group)]
			_groups.extend(addGroups)
			if _groups and _groups == addGroups:
				# track Model deletions
				def checkGroups(tname, myData, tdata):
					if tdata.deleted:
						newGroups = []
						for group in _groups:
							group.update()
							if group.models:
								newGroups.append(group)
						_groups[:] = newGroups
						if not _groups:
							from chimera.triggerSet import ONESHOT
							return ONESHOT
				chimera.triggers.addHandler('Model', checkGroups, None)
		if removeGroups:
			for rg in removeGroups:
				_groups.remove(rg)
	return newGroup

def getGroupOf(model):
	if not _mp:
		raise RuntimeError("Model Panel not yet created")
	from Group import Group
	for item in _mp.items:
		if isinstance(item, Group):
			if model in item.models:
				return item
		elif item == model:
			return None
	raise ValueError("Model not in model panel at all")

def _saveGroups(trigName, myData, sessionFile):
	from Group import Group
	if _mp:
		groups = [g for g in _mp.items if isinstance(g, Group)]
	else:
		groups = _groups
	if not groups:
		return
	from SimpleSession import sessionID
	def _repr(grp):
		strings = []
		from ModelPanel.Group import Group
		for c in grp.components:
			if isinstance(c, Group):
				strings.append("groupCmd(%s)" % _repr(c))
			else:
				try:
					minfo = sessionID(c)
				except:
					minfo = (c.id, c.subid, c.__class__.__name__)
				strings.append("_mpGetModel(%s)" % repr(minfo))
		return "[%s], name=%s" % (", ".join(strings), repr(grp.name))
	print>>sessionFile, """
try:
	def _mpAfterModels():
		def _mpGetModel(info):
			from SimpleSession import modelMap, idLookup
			if isinstance(info, tuple) and len(info) == 3:
				id, subid, className = info
				return [m for m in modelMap[(id, subid)]
					if m.__class__.__name__ == className][0]
			return idLookup(info)
		from ModelPanel import groupCmd
"""
	for grp in groups:
		print>>sessionFile, "\t\tgroupCmd(%s)" % _repr(grp)
	print>> sessionFile, """
	registerAfterModelsCB(_mpAfterModels)
	del _mpAfterModels
except:
	reportRestoreError("Error restoring model panel groups")
"""
from SimpleSession import SAVE_SESSION
chimera.triggers.addHandler(SAVE_SESSION, _saveGroups, None)

from chimera.baseDialog import ModalDialog
class GroupNameDialog(ModalDialog):
	buttons = ("OK",)
	default = "OK"
	#help = 

	def __init__(self, defName):
		self.defName = defName
		ModalDialog.__init__(self)

	def fillInUI(self, parent):
		import Pmw
		self.nameEntry = Pmw.EntryField(parent, labelpos='w', label_text=
			"Group name:")
		self.nameEntry.setentry(self.defName)
		self.nameEntry.component('entry').select_range(0, 'end')
		self.nameEntry.component('entry').focus()
		self.nameEntry.grid(sticky="ew")

	def OK(self):
		self.Cancel(self.nameEntry.getvalue())

def getGroupName(items):
	defName = defaultGroupName(items)
	if chimera.nogui:
		return defName
	master = _mp and _mp.uiMaster().winfo_toplevel() or chimera.tkgui.app
	return GroupNameDialog(defName).run(master)

groupCounter = 0
def defaultGroupName(items):
	global groupCounter
	groupCounter += 1
	name = "Group %d" % groupCounter
	names = set([item.name for item in items])
	if len(names) == 1:
		name = names.pop()
		if not (len(name) == 4 and name[0].isdigit()
				and name[1:].isalnum()):
			# not PDB ID
			name = plural(name)
	else:
		import os.path
		prefix = os.path.commonprefix(names)
		if len(prefix) > 2:
			for n in names:
				if len(n) > len(prefix) and n[len(prefix)].isalpha():
					break
			else:
				name = plural(prefix)
	return name

def plural(text):
	if not text:
		return text
	if text[-1] == 's':
		if text[-2:] == "ss":
			return text + "es"
		else:
			return text
	if text[-1] == 'h':
		if text[-2:-1] in "cs":
			return text + "es"
		else:
			return text + "s"
	if text[-1] in "oxz":
		return text + "es"
	return text + "s"

def noteCmd(models):
	from noteDialog import NoteDialog
	NoteDialog(models)

def renameCmd(items):
	from renameDialog import RenameDialog
	RenameDialog(items)

def viewerCmd(models):
	viewer.viewAll()

def selectCmd(models):
	sel = chimera.selection.ItemizedSelection()
	sel.add(models)
	chimera.tkgui.selectionOperation(sel)
	
	from chimera import openModels, viewer, update
	shown = {}
	for m in openModels.list():
		shown[m] = m.display
		if m in models:
			m.display = 1
		else:
			m.display = 0
	update.checkForChanges()
	viewer.viewAll()
	if chimera.openModels.cofrMethod != chimera.OpenModels.Independent:
		openModels.cofrMethod = openModels.CenterOfView
		viewer.clipping = True

	for m,disp in shown.items():
		m.display = disp
	update.checkForChanges()
	
def showAllAtomsCmd(models):
	for m in models:
		if not hasattr(m, 'atoms') or not hasattr(m, 'bonds'):
			continue
		m.display = 1
		for a in m.atoms:
			a.display = 1

def surfCmd(models, category):
	import Midas
	mols = filter(lambda m: isinstance(m, chimera.Molecule), models)
	Midas.surfaceNew(category, models=mols)
	for m in mols:
		for a in m.atoms:
			if a.surfaceCategory == category:
				a.surfaceDisplay = 1

def _getModels(item):
	# VolumeModel has a 'models' attr(!), so...
	if isinstance(item, Group):
		return item.models
	return [item]

import threading

class DistanceCompareThread(threading.Thread):
	def __init__(self, tracerPath, path, mrcPath, skeletonPath, threshold, pdbPath, outputQueue, outPath):
		self.tracerPath = tracerPath
		self.path = path
		self.pdbPath = pdbPath
		self.mrcPath = mrcPath
		self.skeletonPath = skeletonPath
		self.threshold = threshold
		self.outputQueue = outputQueue
		self.outPath = outPath
		threading.Thread.__init__(self)

	def run(self):
		if (self.skeletonPath == ""):
			self.skeletonPath = "Empty"
		if (self.threshold == ""):
			self.threshold = "Empty"
		if (self.outPath == ""):
			self.outPath = "Empty"
		if platform.system() == 'Windows':
			arguments = '"' + self.mrcPath +  '" "'  + self.skeletonPath + '" "' + self.threshold + '" "' + self.outPath + '"'
		else:
			arguments = '"' + self.mrcPath +  '" "'  + self.skeletonPath + '" "' + self.threshold + '" "' + self.outPath + '"'
		# for i in self.chains:
		# 	arguments += " " + i
		#self.outputBox.insert(Tkinter.END, "Running: " + self.cutDensityPath + "\n with arguments: " + arguments + "\n\n")
		#self.outputQueue.put("Running: " + self.tracerPath + "\n" + arguments + "\n\n")
		#process = subprocess.Popen('exec "' + self.tracerPath + " " + arguments + '"', stdout=subprocess.PIPE, shell=True)
		process = subprocess.Popen(self.tracerPath + " " + arguments, stdout=subprocess.PIPE, shell=True)
		#process = subprocess.Popen([self.tracerPath, self.pdbPath, pdbID, stdout=subprocess.PIPE, shell=False)

		while True:
			line = process.stdout.readline()

			if line != '':
				self.outputQueue.put(line.rstrip() + "\n")
				#self.outputBox.insert(Tkinter.END, line.rstrip() + "\n")
			else:
				break

		self.outputQueue.put("Process finished\n")

class TracerThread(threading.Thread):
	def __init__(self, tracerPath, path, mrcPath, skeletonPath, threshold, analysis, pdbPath, openCut, outputQueue):
		self.openCut = openCut
		self.tracerPath = tracerPath
		self.path = path
		self.pdbPath = pdbPath
		self.mrcPath = mrcPath
		self.skeletonPath = skeletonPath
		self.threshold = threshold
		self.outputQueue = outputQueue
		self.analysis = analysis
		threading.Thread.__init__(self)

	def run(self):
		if platform.system() == 'Windows':
			arguments = '"' + self.path + '"\ "' + self.mrcPath + '" "' + self.skeletonPath + '" ' + str(self.threshold) + ' ' + str(self.analysis) + ' "' + self.pdbPath + '"'
		else:
			arguments = '"' + self.path + '" "' + self.mrcPath + '" "' + self.skeletonPath + '" ' + str(self.threshold) + ' ' + str(self.analysis) + ' "' + self.pdbPath + '"'
		# for i in self.chains:
		# 	arguments += " " + i

		#self.outputBox.insert(Tkinter.END, "Running: " + self.cutDensityPath + "\n with arguments: " + arguments + "\n\n")
		self.outputQueue.put("Running: " + self.tracerPath + "\n" + arguments + "\n\n")
		#process = subprocess.Popen('exec "' + self.tracerPath + " " + arguments + '"', stdout=subprocess.PIPE, shell=True)
		process = subprocess.Popen(self.tracerPath + " " + arguments, stdout=subprocess.PIPE, shell=True)
		#process = subprocess.Popen([self.tracerPath, self.pdbPath, pdbID, stdout=subprocess.PIPE, shell=False)

		while True:
			line = process.stdout.readline()

			if line != '':
				self.outputQueue.put(line.rstrip() + "\n")
				#self.outputBox.insert(Tkinter.END, line.rstrip() + "\n")
			else:
				break

		self.outputQueue.put("Process finished\n")

		if self.openCut:
			import glob
			import time

			#1FLP_thr_0.4_outFiles

			pdbID = os.path.basename(self.mrcPath)

			outputPath = self.mrcPath + "_thr_" + self.threshold + "_outFiles"
			#self.outputQueue.put(outputPath + "\n")

			openModelsBefore = list(chimera.openModels.list())

			#os.chdir(outputPath)
			for file in glob.glob(outputPath + os.sep + '*.pdb'):
				self.outputQueue.put("Opening: " + file + "\n")
				chimera.openModels.open(file)
				time.sleep(0.5)

			#1FLP_SHTestimate_final.mrc

			#currently /home/taylor/Documents/Research/1FLP/1FLP_thr_0.41_outFiles/1FLP_SHTestimate_final.mrc

			#group all as self.pdbID + "_thr_" + self.threshold
			#get list of open items before and after opening output files, only group the new ones
			finalFile = outputPath + os.sep + pdbID + "_SHTestimate_final.mrc"
			chimera.openModels.open(finalFile)

			openModelsNow = chimera.openModels.list()

			outputModels = set(openModelsNow) - set(openModelsBefore)

			groupName = pdbID + "_thr_" + self.threshold
			ModelPanel.groupCmd(outputModels, groupName)

			#self.outputQueue.put(finalFile + "\n")

			self.outputQueue.put("Opened output files\n")

class TwisterThread(threading.Thread):
	def __init__(self, twisterPath, path, mrcPath, skeletonPath, threshold, analysis, pdbPath, openCut, outputQueue):
		self.openCut = openCut
		self.twisterPath = twisterPath
		self.path = path
		self.pdbPath = pdbPath
		self.mrcPath = mrcPath
		self.skeletonPath = skeletonPath
		self.threshold = threshold
		self.outputQueue = outputQueue
		self.analysis = analysis
		threading.Thread.__init__(self)

	def run(self):
		if platform.system() == 'Windows':
			arguments = '"' + self.pdbPath + '" "' + self.mrcPath + '" ' + str(self.threshold);
		else:
			arguments = '"' + self.pdbPath + '" "' + self.mrcPath + '" ' + str(self.threshold);
		# for i in self.chains:
		# 	arguments += " " + i

		#self.outputBox.insert(Tkinter.END, "Running: " + self.cutDensityPath + "\n with arguments: " + arguments + "\n\n")
		self.outputQueue.put("Running: " + self.twisterPath + "\n" + arguments + "\n\n")
		#process = subprocess.Popen('exec "' + self.tracerPath + " " + arguments + '"', stdout=subprocess.PIPE, shell=True)
		process = subprocess.Popen(self.twisterPath + " " + arguments, stdout=subprocess.PIPE, shell=True)
		#process = subprocess.Popen([self.tracerPath, self.pdbPath, pdbID, stdout=subprocess.PIPE, shell=False)

		while True:
			line = process.stdout.readline()

			if line != '':
				self.outputQueue.put(line.rstrip() + "\n")
				#self.outputBox.insert(Tkinter.END, line.rstrip() + "\n")
			else:
				break

		self.outputQueue.put("Process finished\n")

		if self.openCut:
			import glob
			import time

			#1FLP_thr_0.4_outFiles

			pdbID = os.path.basename(self.pdbPath)

			outputPath = self.mrcPath + "twist" #"_thr_" + self.threshold + "_outFiles"
			self.outputQueue.put("Output files written to " + outputPath +" \n")

			openModelsBefore = list(chimera.openModels.list())

			#4ZQQsht37_trans_2_orient_10.pdb
			for file in glob.glob(outputPath + os.sep + pdbID + '*[0-9].pdb'): #_trans_[0-9]_orient
				if not "CAchain" in file:
					self.outputQueue.put("Opening: " + file + "\n")
					chimera.openModels.open(file)
					time.sleep(0.5)

			#1FLP_SHTestimate_final.mrc

			#currently /home/taylor/Documents/Research/1FLP/1FLP_thr_0.41_outFiles/1FLP_SHTestimate_final.mrc

			#group all as self.pdbID + "_thr_" + self.threshold
			#get list of open items before and after opening output files, only group the new ones

			openModelsNow = chimera.openModels.list()

			outputModels = set(openModelsNow) - set(openModelsBefore)

			groupName = pdbID + "sht_thr_" + self.threshold
			ModelPanel.groupCmd(outputModels, groupName)

			#self.outputQueue.put(finalFile + "\n")

			self.outputQueue.put("Opened output files\n")
						
		
			
			
#dialogs.register(SSETracerDialog.name, SSETracerDialog)