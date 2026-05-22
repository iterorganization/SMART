#!/usr/bin/env python3
"""
Wrapper around iWrap 1.0.0 that adds an `imas.IDSName` shim for
compatibility with IMAS-Python 2.x.

iWrap 1.0.0 (the version on SDCC at the time of writing) calls
`list(imas.IDSName)` inside `iwrap/common/utils.py:get_all_ids_names()`
as a last-resort fallback when the `data_dictionary` Python package is
absent.  IMAS-Python 2.x does not expose an `IDSName` enum; the
replacement is `imas.IDSFactory().ids_names()`.  We synthesise an enum
on the fly so iWrap's existing code path keeps working.

Remove this wrapper (and the Makefile reference to it) once iWrap is
updated to support IMAS-Python 2.x natively.
"""
import sys
import imas

if not hasattr(imas, "IDSName"):
    from enum import Enum
    imas.IDSName = Enum(
        "IDSName",
        {name: name for name in imas.IDSFactory().ids_names()},
    )

# Match iWrap's own bin/iwrap entry script
from bin.run import cmd_line
sys.exit(cmd_line())
