#!/usr/bin/env python
# Copyright 2025 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Utilities for processing command line arguments.
"""
from __future__ import division, print_function
import sys, itertools

def _get_arg_pos(string, silent = False) :
    """
    Return the position of string in sys.argv, or -1 if not present.
    If not silent, print a message that the argument is being used.
    """
    for ii, arg in enumerate(sys.argv) :
        if string != arg :
            continue
        if not silent :
            print('Using command line argument: ' + string, file = sys.stderr)
        return ii
    return -1

def check_arg(string, remove = False, silent = False) :
    """
    Return true if string is one of the arguments.
    Remove it from sys.argv if requested.
    If not silent, print a message that the argument is being used.
    """
    pos = _get_arg_pos(string, silent)
    if pos < 0 :
        return False
    if remove :
        del sys.argv[pos]
    return True

def get_associated_arg(string, default = '', remove = False, silent = False,
                       numArgs = 1, mappers = None) :
    """
    Look for an argument named string, and return the subsequent argument(s) as
        a single value if numArgs == 1, otherwise a list.
    If string isn't there, return default
    If string is there but is last, alert the user and stop.
    If remove is true, remove the string and associated argument from sys.argv
    If mappers is not None, it is a function that maps arguments (e.g., int) or a
        sequence of such functions, which will be applied to (non-default) output(s).
    If not silent, print a message that the argument is being used.
    """
    pos = _get_arg_pos(string, silent)
    if pos < 0 :
        return default
    if pos >= len(sys.argv) - numArgs :
        if numArgs == 1 :
            print('Argument %s requires a value.' % string, file = sys.stderr)
        else :
            print('Argument %s requires %d values.' % (string, numArgs), file = sys.stderr)
        raise SystemExit(1)
    if not silent :
        if numArgs == 1 :
            print('   with value %s' % sys.argv[pos + 1], file = sys.stderr)
        else :
            print('   with values %s' % sys.argv[pos + 1 : pos + numArgs + 1],
                  file = sys.stderr)
    assocArgs = list(sys.argv[pos + 1 :  pos + numArgs + 1])
    if remove :
        del sys.argv[pos : pos + numArgs + 1]
    if mappers is not None :
        # If mappers is not iterable, replace it with [mappers]
        try :
            iter(mappers)
        except TypeError :
            mappers = [mappers]

        for argInd, mapper in zip(range(len(assocArgs)), itertools.cycle(mappers)) :
            if mapper is not None :
                assocArgs[argInd] = mapper(assocArgs[argInd])
    return assocArgs[0] if numArgs == 1 else assocArgs

get_associated_arg_str = get_associated_arg

def get_associated_arg_int(*pArgs, **kArgs) :
    kArgs['mappers'] = int
    return get_associated_arg(*pArgs, **kArgs)

def get_associated_arg_float(*pArgs, **kArgs) :
    kArgs['mappers'] = float
    return get_associated_arg(*pArgs, **kArgs)

def assert_num_args(numArgsRequired, correctUsageString, exact = False,
                    maxNumAllowed = None, silent = False, forceFailure = False) :
    """
    Check that number of arguments is at least numArgsRequired and if not then print
        usage string (if not silent) and stop.
    If maxNumAllowed is not None, require that there are at most that many arguments.
    If exact is True, require that there are exactly numArgsRequired arguments (equivalent
        to setting maxNumAllowed to numArgsRequired).
    If forceFailure, always print usage string and stop even if there are the correct
        number of arguments. This is useful if something else is wrong with arguments.
    """
    numArgs = len(sys.argv) - 1
    if (numArgs < numArgsRequired or
            (exact and numArgs > numArgsRequired) or
            (maxNumAllowed is not None and numArgs > maxNumAllowed) or
            forceFailure) :
        if not silent :
            print('Usage: ' + correctUsageString, file = sys.stderr)
        raise SystemExit(1)
