#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re

def main():
    if(len(sys.argv) < 2):
        print("Usage: gfanalyse_gfortran_ast.py <ast-file-1> ... <ast-file-N>")
        print("       For generating the abstract syntax tree (ast) run gfortran")
        print('       with "-fdump-fortran-original" and redirect output to file.')
        print('       This can be achieved by putting "FC_SAVELOG = TRUE" in the cp2k arch-file.')
        sys.exit(1)

    log_files = sys.argv[1:]

    for fn in log_files:
        process_log_file(fn)

#===============================================================================
def process_log_file(fn):
    log = open(fn).read()
    lines = log.split("\n")
    curr_procedure = None
    curr_symbol = None

    re_procedure = re.compile(r"^\s+procedure name = (.*)$")
    re_symbol    = re.compile(r"^\s+symtree.* symbol: '([^']+)'.*$")
    re_attr      = re.compile(r"^\s+attributes: (.*)$")
    re_open      = re.compile(r"^\s+OPEN")
    re_close     = re.compile(r"^\s+CLOSE")

    open_close_exceptions = ["cp_files", "machine_gfortran", "machine",]

    for line in lines:
        m = re_procedure.match(line)
        if(m):
            curr_procedure = m.group(1)
            continue

        m = re_symbol.match(line)
        if(m):
            curr_symbol = m.group(1)
            continue

        m = re_attr.match(line)
        if(m and curr_procedure and ("IMPLICIT-SAVE" in line) and ("PARAMETER" not in line)):
            msg  = fn+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-SAVE'
            print(msg)

        m = re_attr.match(line)
        if(m and curr_procedure and ("IMPLICIT-TYPE" in line)):
            msg  = fn+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-TYPE'
            print(msg)

        m = re_open.match(line)
        if(m and curr_procedure and (fn.rsplit(".",1)[0] not in open_close_exceptions)):
            msg = fn+': direct call to OPEN in procedure "'+curr_procedure+"'"
            print(msg)

        m = re_close.match(line)
        if(m and curr_procedure and (fn.rsplit(".",1)[0] not in open_close_exceptions)):
            msg = fn+': direct call to CLOSE in procedure "'+curr_procedure+"'"
            print(msg)

#===============================================================================
main()

#EOF