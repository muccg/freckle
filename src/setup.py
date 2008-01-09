#!/usr/bin/env python

from distutils.core import setup

setup(name='freckle',
	version='0.1',
	description='Freckle dotplot utility',
	author='Crispin Wellington',
	author_email='cwellington@ccg.murdoch.edu.au',
	packages=['pyfreckle'],
	scripts=['freckle','DotPlot.py'],
	)