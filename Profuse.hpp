/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      Profuse.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Global #defines for the Profuse library.
##
#******************************************************************************
#*
#*    This file is part of profuse, a suite of programs for working with
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (LGPL v3), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
#*    profuse is free software: you can redistribute it and/or modify it under
#*    the terms of the GNU Lesser Public License as published by the Free
#*    Software Foundation, either version 3 of the License, or (at your option)
#*    any later version.
#*
#*    profuse is distributed in the hope that it will be useful, but WITHOUT
#*    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#*    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for
#*    more details.
#*
#*    You should have received a copy of the GNU Lesser Public License along
#*    with profuse.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROFUSE_HPP__
#define __GALOSH_PROFUSE_HPP__

#include "Prolific.hpp"

// ProfileTrainer
#define GALOSH_INLINE_TRAIN
#define GALOSH_INLINE_TRAIN_PERCENT_CHANGE inline

// ProfileGibbs
#define GALOSH_INLINE_GIBBS

// SimulationStudy
#define GALOSH_INLINE_SIMULATIONSTUDY_START

#endif // __GALOSH_PROFUSE_HPP__
