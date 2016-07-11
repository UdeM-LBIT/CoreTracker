# CoreTracker Copyright (C) 2016  Emmanuel Noutahi
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import warnings
from corefile import CoreFile
from collections import defaultdict
from functools import partial

warnings.filterwarnings("ignore")

import utils
from utils import SequenceLoader, SequenceSet, ReaGenomeFinder

def makehash(depth, type):
    if depth == 0:
        return defaultdict(type)
    else:
        return defaultdict(partial(makehash, depth-1, type))

__all__ = ['utils', 'SequenceLoader', 'SequenceSet', 'CoreFile',
            'ReaGenomeFinder', 'makehash']
