# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: writePDBdialog.py 38071 2013-01-03 23:44:32Z pett $

"""dialog for writing models as PDB files"""

from OpenSave import SaveModeless
from chimera import replyobj, preferences, triggers, openModels, Molecule

class WritePDBdialog(SaveModeless):
	keepShown = SaveModeless.default
	help = "UsersGuide/savemodel.html"
	name = "write PDB"

	def __init__(self):
		self.prefs = preferences.addCategory("write PDB dialog",
				preferences.HiddenCategory,
				optDict={"multiSavePDB": "multiple"})
		self.haveTraj = False
		SaveModeless.__init__(self, clientPos='s', clientSticky='ewns',
			filters=[("PDB", "*.pdb", ".pdb")])
		openModels.addAddHandler(self._modelsChange, None),
		openModels.addRemoveHandler(self._modelsChange, None)
		self._modelsChange()

	def configure(self, models=None, refreshList=True, selOnly=None):
		if models is not None:
			if len(models) > 1:
				name = "Multiple Models "
			elif models:
				name = models[0].name + " "
			else:
				name = ""
			self._toplevel.title("Save %sas PDB File" % name)
			if refreshList:
				self.modelList.setvalue(models)

			self._trajCheck()

			if len(models) > 1:
				self.multiSaveMenu.grid(row=self._msmRow,
						column=0, sticky='w')
			else:
				self.multiSaveMenu.grid_forget()
		if selOnly is not None:
			self.selOnlyVar.set(selOnly)

	def fillInUI(self, parent):
		SaveModeless.fillInUI(self, parent)
		row = 0

		from chimera.widgets import MoleculeScrolledListBox, \
							ModelOptionMenu
		self.modelList = MoleculeScrolledListBox(self.clientArea,
			labelpos='w', label_text="Save models:",
			listbox_selectmode='extended',
			selectioncommand=lambda: self.configure(
			self.modelList.getvalue(), refreshList=False))
		self.modelList.grid(row=row, column=0, sticky='nsew')
		self.clientArea.rowconfigure(row, weight=1)
		self.clientArea.columnconfigure(0, weight=1)
		row += 1

		import Tkinter, Pmw
		self.dispOnlyVar = Tkinter.IntVar(parent)
		self.dispOnlyVar.set(False)
		Tkinter.Checkbutton(self.clientArea, variable=self.dispOnlyVar,
			text="Save displayed atoms only").grid(row=row,
			column=0, sticky='w')
		row += 1

		self.selOnlyVar = Tkinter.IntVar(parent)
		self.selOnlyVar.set(False)
		Tkinter.Checkbutton(self.clientArea, variable=self.selOnlyVar,
			text="Save selected atoms only").grid(row=row,
			column=0, sticky='w')
		row += 1

		self.saveRelativeVar = Tkinter.IntVar(parent)
		self.saveRelativeVar.set(False)
		self.relativeFrame = f = Tkinter.Frame(self.clientArea)
		Tkinter.Checkbutton(f, variable=self.saveRelativeVar,
			text="Save relative to model:"
			).grid(row=0, column=0, sticky='e')
		self.relModelMenu = ModelOptionMenu(f)
		self.relModelMenu.grid(row=0, column=1, sticky='w')
		self.saveUntransformedVar = Tkinter.IntVar(parent)
		self.saveUntransformedVar.set(True)
		self.untransformedButton = Tkinter.Checkbutton(self.clientArea,
					variable=self.saveUntransformedVar,
					text="Use untransformed coordinates")
		self._rfRow = row
		row += 1

		self.frameSave = Pmw.OptionMenu(self.clientArea,
			labelpos='w', label_text="Save",
			initialitem="current frame",
			items=["current frame", "all frames"])
		# not always shown; remember row number
		self._fsRow = row
		row += 1
			
		from chimera import dialogs
		self.labelMap = {
			"single": "a single file",
			"multiple":
				"multiple files [file name must contain $name or $number]"
		}
		preferred = self.labelMap[self.prefs["multiSavePDB"]]
		self.multiSaveMenu = Pmw.OptionMenu(self.clientArea,
			labelpos='w', label_text="Save multiple models in",
			initialitem=preferred, items=self.labelMap.values())
		# not always shown; remember row number
		self._msmRow = row
		row += 1

	def map(self, *args):
		self.handler = triggers.addHandler("CoordSet",
							self._trajCheck, None)
		self._trajCheck()

	def unmap(self, *args):
		triggers.deleteHandler("CoordSet", self.handler)
		self.handler = None

	def Apply(self):
		from chimera import dialogs

		selOnly = self.selOnlyVar.get()

		paths = self.getPaths()
		if not paths:
			replyobj.error('No save location chosen.\n')
			return
		path = paths[0]
		models = self.modelList.getvalue()
		if not models:
			replyobj.error("No models chosen to save.\n")
			return
		if len(openModels.listIds()) > 1:
			if self.saveRelativeVar.get():
				relModel = self.relModelMenu.getvalue()
			else:
				relModel = None
		else:
			if self.saveUntransformedVar.get():
				relModel = models[0]
			else:
				relModel = None
		if self.haveTraj and self.frameSave.getvalue() == "all frames":
			allFrames=True
		else:
			allFrames=False
			
		import Midas
		import os
		# path will be encoded if the OS doesn't support Unicode
		# file names, so decode before printing
		if os.path.supports_unicode_filenames:
			printablePath = path
		else:
			printablePath = path.decode('utf8')
		if len(models) < 2:
			replyobj.status("Writing %s to %s" %
							(models[0].name, printablePath))
			Midas.write(models, relModel, path,
					allFrames=allFrames,
					dispOnly=self.dispOnlyVar.get(),
					selOnly=selOnly)
			replyobj.status("Wrote %s to %s" %
							(models[0].name, printablePath))
			return

		saveOpt = self.multiSaveMenu.getvalue()
		for key, value in self.labelMap.items():
			if saveOpt == value:
				break
		self.prefs["multiSavePDB"] = key

		# write multiple models to multiple files
		if key == "multiple":
			if "$name" not in path and "$number" not in path:
				self.enter()
				from chimera import UserError
				raise UserError("You must provide a file name that contains"
					" either $name or $number, which will be changed to"
					" the model name or model number (respectively) in the"
					" final saved file name")
			for m in models:
				modelPath = path.replace("$name", m.name).replace(
					"$number", m.oslIdent()[1:])
				if os.path.supports_unicode_filenames:
					printablePath = modelPath
				else:
					printablePath = modelPath.decode('utf8')
				replyobj.status("Writing %s (%s) to %s" %
					(m.name, m.oslIdent(), printablePath))
				# since we can write a _lot_ of files here,
				# use temporary=True to avoid washing out
				# the Rapid Access list
				Midas.write(m, relModel, modelPath,
					allFrames=allFrames,
					dispOnly=self.dispOnlyVar.get(),
					selOnly=selOnly, temporary=True)
				replyobj.status("Wrote %s (%s) to %s" %
					(m.name, m.oslIdent(), printablePath))
			return

		# write multiple models to single file
		replyobj.status("Writing multiple models to %s" % printablePath)
		Midas.write(models, relModel, path, allFrames=allFrames,
					dispOnly=self.dispOnlyVar.get(),
					selOnly=selOnly)
		replyobj.status("Wrote multiple models to %s" % printablePath)

	def _modelsChange(self, *args):
		# can't query listbox, since it hangs off of same trigger
		if len(openModels.listIds()) > 1:
			self.untransformedButton.grid_forget()
			self.relativeFrame.grid(row=self._rfRow,
						column=0, sticky='w')
		else:
			self.relativeFrame.grid_forget()
			self.untransformedButton.grid(row=self._rfRow,
						column=0, sticky='w')
	def _trajCheck(self, *args):
		haveTraj = False
		for m in self.modelList.getvalue():
			if len(m.coordSets) > 1:
				haveTraj = True
				break
		if haveTraj == self.haveTraj:
			return
		self.haveTraj = haveTraj
		if self.haveTraj:
			self.frameSave.grid(row=self._fsRow, column=0,
								sticky='w')
		else:
			self.frameSave.grid_forget()

from chimera import dialogs
dialogs.register(WritePDBdialog.name, WritePDBdialog)
