# --- UCSF Chimera Copyright ---
# Copyright (c) 2000-2006 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: noteDialog.py 26655 2009-01-07 22:02:30Z gregc $

import chimera
from chimera.baseDialog import ModelessDialog

class RenameDialog(ModelessDialog):
	oneshot = True
	title = "Rename"
	buttons = ("OK", "Cancel")
	default = "OK"
	help="UsersGuide/modelpanel.html#rename"

	def __init__(self, items):
		self.items = items
		ModelessDialog.__init__(self)

	def fillInUI(self, parent):
		import Tkinter
		from Group import Group
		contents = None
		for item in self.items:
			if isinstance(item, Group):
				if contents is None:
					contents = "groups"
				elif contents != "groups":
					contents = "models/groups"
					break
			else:
				if contents is None:
					contents = "models"
				elif contents != "models":
					contents = "models/groups"
					break
		names = set([item.name for item in self.items])
		if len(names) == 1:
			mname = names.pop()
		else:
			mname = "new name"
		Tkinter.Label(parent, text="Rename to:").grid(row=0, sticky='w')
		self.nameVar = Tkinter.StringVar(parent)
		self.nameVar.set(mname)
		entry = Tkinter.Entry(parent, exportselection=False,
						textvariable=self.nameVar)
		entry.focus_set()
		entry.selection_range(0, Tkinter.END)
		entry.icursor(Tkinter.END)
		entry.grid(row=1, sticky='ew')

		self.renameModelsVar = Tkinter.IntVar(parent)
		self.renameModelsVar.set(contents != "groups")
		Tkinter.Checkbutton(parent, text="Rename models",
			variable=self.renameModelsVar).grid(row=2, column=0, sticky="w")

		self.renameGroupsVar = Tkinter.IntVar(parent)
		self.renameGroupsVar.set(contents == "groups")
		ck = Tkinter.Checkbutton(parent, text="Rename groups",
			variable=self.renameGroupsVar)
		ck.grid(row=3, column=0, sticky="w")
		if contents == "models":
			ck.configure(state="disabled")

	def Apply(self):
		mname = self.nameVar.get()
		from Group import Group
		needsPoking = False
		if self.renameGroupsVar.get():
			for item in self.items:
				if isinstance(item, Group):
					item.name = mname
					needsPoking = True
		if self.renameModelsVar.get():
			for item in self.items:
				if item.__destroyed__:
					continue
				if isinstance(item, Group):
					for m in item.models:
						m.name = mname
				else:
					item.name = mname
				needsPoking = False
		if needsPoking:
			from base import _mp
			_mp._fillTable()
