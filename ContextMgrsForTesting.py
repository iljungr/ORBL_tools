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
Context managers useful for testing, such as temporarily changing and restoring os.environ
    or sys.argv, or temporarily capturing stdout or stderr.
"""
from __future__ import division, print_function
import sys, os, copy
if sys.version_info[0] < 3 : # pragma: no cover
    from StringIO import StringIO
else : # pragma: no cover
    from io import StringIO

class StreamCatcher(object) :
    """
    Context manager that grabs all output to sys.stdout, sys.stderr, or both and puts
        it (them) in a string (or strings). The string(s) can be accessed even after the
        context manager has exited.
    Usage:
        with StreamCatcher(kind) as streamCatcher : # kind in ['out', 'err', 'both']
            do stuff
            check that streamCatcher.buffer('out') or streamCatcher.buffer('err')
                contains what is expected.
            Optionally, streamCatcher.clear('out') or streamCatcher.clear('err') and then
                do more stuff and more checking
        do more stuff with streamCatcher.buffer('out') or streamCatcher.buffer('err')
    If kind is not 'both', buffer() may be used instead of buffer('out') or buffer('err'),
        and clear() may be used insteaad of clear('out') or clear('err').
    Implementation note: to simplify the code, we always save both stderr and stdout,
        but we only suppress output to the two streams if asked to do so.
    """
    def __init__(self, kind = 'both') :
        assert kind in ['out', 'err', 'both']
        self._kind = kind # Please don't change after this initialization!
    def __enter__(self) :
        self.savedStreams = [sys.stdout, sys.stderr]
        self._buffers     = [StringIO(), StringIO()] # stdout, stderr
        self.savedBuffers = None
        if self._kind != 'err' :
            sys.stdout = self._buffers[0]
        if self._kind != 'out' :
            sys.stderr = self._buffers[1]
        return self
    def __exit__(self, exc_type, exc_val, exc_tb) :
        if self._kind != 'err' :
            sys.stdout = self.savedStreams[0]
        if self._kind != 'out' :
            sys.stderr = self.savedStreams[1]
        self.savedBuffers = [self._buffers[0].getvalue(), self._buffers[1].getvalue()]
        for _buf in self._buffers :
            _buf.close()
    def _check_which(self, which) :
        """Check that 'which' argument is consistent with self._kind """
        assert which in ['out', 'err', 'both', 'whichever']
        assert which != 'whichever' or self._kind != 'both'
        assert which != 'err' or self._kind != 'out'
        assert which != 'out' or self._kind != 'err'
        assert which != 'both' or self._kind == 'both'
    def buffer(self, which = 'whichever') : # which in 'out', 'err', 'both', 'whichever'
        """
        Return output to stdout or stderr since self was created or last cleared.
        If which is 'whichever', use the one specified on creation, which may not be both.
        If which is 'both', return a pair of out and err
        """
        self._check_which(which)
        if which == 'whichever' :
            which = self._kind
        def get_buffer(ind) :
            return (self._buffers[ind].getvalue() if self.savedBuffers is None else
                    self.savedBuffers[ind])
        return (get_buffer(1) if which == 'err' else get_buffer(0) if which == 'out' else
                (get_buffer(0), get_buffer(1)))
    def clear(self, which = 'whichever') : # which in 'out', 'err', 'both', 'whichever'
        """Clear the buffer so that future writes will be on a clean slate."""
        self._check_which(which)
        if which == 'whichever' :
            which = self._kind
        for ind in [0] if self._kind == 'out' else [1] if self._kind == 'err' else [0, 1]:
            if self.savedBuffers is None :
                self._buffers[ind].truncate(0)
                self._buffers[ind].seek(0)
            else : # Not sure why you'd even clear after exit, but if you do this works
                self.savedBuffers[ind] = ''

class StderrCatcher(StreamCatcher) :
    """
    Usage:
        with StderrCatcher() as errCatcher :
            do stuff
            check that errCatcher.buffer() contains what is expected.
            Optionally, errCatcher.clear() and then  do more stuff and more checking
        do more stuff with errCatcher.buffer()
    """
    def __init__(self) :
        StreamCatcher.__init__(self, 'err')

class StdoutCatcher(StreamCatcher) :
    """
    Usage:
        with StderrCatcher() as errCatcher :
            do stuff
            check that errCatcher.buffer() contains what is expected.
            Optionally, errCatcher.clear() and then  do more stuff and more checking
        do more stuff with errCatcher.buffer()
    """
    def __init__(self) :
        StreamCatcher.__init__(self, 'out')

class EnvSaver(object) :
    """
    Context manager for temporary changes to os.environ that guarantees it is restored.
    Usage:
        with EnvSaver() :
            do stuff that could change environment
        # Now environment has been restored.
    """
    def __enter__(self) :
        self.backupEnviron = copy.deepcopy(os.environ)
        return self
    def __exit__(self, exc_type, exc_val, exc_tb) :
        os.environ.clear()
        os.environ.update(self.backupEnviron)

class SysArgSaver(object) :
    """
    Context manager for temporary changes to sys.argv that guarantees it is restored.
    Usage:
        with SysArgSaver() :
            do stuff that could change sys.argv
        # Now sys.argv has been restored.
    """
    def __enter__(self) :
        self.backupSysArgv = copy.deepcopy(sys.argv)
        return self
    def __exit__(self, exc_type, exc_val, exc_tb) :
        sys.argv = self.backupSysArgv