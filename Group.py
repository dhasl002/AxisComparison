# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: __init__.py 26655 2009-01-07 22:02:30Z gregc $

groupCounter = 0

class Group:
	def __init__(self, components, name):
		self.components = components
		self.name = name

	def getModels(self):
		models = []
		for c in self.components:
			if isinstance(c, Group):
				models.extend(c.models)
			else:
				models.append(c)
		return models
	models = property(getModels)

	def update(self):
		newComponents = []
		for c in self.components:
			if isinstance(c, Group):
				c.update()
				if c.models:
					newComponents.append(c)
			else:
				if not c.__destroyed__:
					newComponents.append(c)
		self.components = newComponents

	localAttrs = ("name", "components")
	def __getattr__(self, name):
		if name in Group.localAttrs:
			return self.__dict__["_" + name]
		vals = set()
		for c in self.components:
			v = getattr(c, name)
			if type(v) == set:
				vals |= v
			else:
				vals.add(v)
		if len(vals) == 1:
			return vals.pop()
		return GroupAttr(vals)

	def __setattr__(self, name, val):
		if name in Group.localAttrs:
			self.__dict__["_" + name] = val
		else:
			for c in self.components:
				setattr(c, name, val)

	def __eq__(self, group):
		# need explicit __eq__ because otherwise a moribund
		# group with no components is problematic
		if not isinstance(group, Group):
			return False
		return self.components == group.components

	def __str__(self):
		vals = set([str(c) for c in self.components])
		if len(vals) == 1:
			return vals.pop()
		return "---multiple---"
	
	def __repr__(self):
		return "%s(%s)" % (self.__class__.__name__, repr(self.components))

class GroupAttr:
	def __init__(self, vals):
		self.vals = vals

	localAttrs = ("vals",)
	def __getattr__(self, name):
		if name in GroupAttr.localAttrs:
			return self.__dict__["_" + name]
		subVals = set()
		for v in self.vals:
			subVals.add(getattr(v, name))
		if len(subVals) == 1:
			return subVals.pop()
		return GroupAttr(subVals)

	def __setattr__(self, name, val):
		if name in GroupAttr.localAttrs:
			self.__dict__["_" + name] = val
		else:
			for v in self.vals:
				setattr(v, name, val)

	def __call__(self, *args, **kw):
		subVals = set([v(*args, **kw) for v in self.vals])
		if len(subVals) == 1:
			return subVals.pop()
		return GroupAttr(subVals)

	def __lt__(self, other):
		return min(self.vals) < other

	def __gt__(self, other):
		return min(self.vals) > other

	def __eq__(self, other):
		if not isinstance(other, GroupAttr):
			return False
		return self.vals == other.vals

	def __ne__(self, other):
		return not self == other

	def __hash__(self):
		return id(self)

	def __str__(self):
		subVals = set([str(v) for v in self.vals])
		if len(subVals) == 1:
			return subVals.pop()
		raise ValueError("str() called on non-homogenous group attr: " +
			", ".join([str(v) for v in self.vals]))

	def __str__(self):
		return self._stringize(str)

	def __unicode__(self):
		return self._stringize(unicode)

	def _stringize(self, converter):
		subVals = set([converter(v) for v in self.vals])
		if len(subVals) == 1:
			return subVals.pop()
		return "%s-%s" % (converter(min(self.vals)), converter(max(self.vals)))
		#return "---multiple---"

	def __repr__(self):
		return "%s(%s)" % (self.__class__.__name__, repr(self.vals))
