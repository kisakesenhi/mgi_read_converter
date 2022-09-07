# mgi_read_converter

A command line tool that converts sequencing read read headers from MGI format into illumina format.

This is applicable for phred+33 scored base qualities, since some software might assume base qualities are scored phred+64.

MGI read name format => @readname/{pair}   -->  @F345455..5654/1

Illumina read format => @readname:0:0:0:0:0:0 {pair}:N:0:1 --> @F345455..5654 1:N:0:1



### GNU GPL v3
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)    
`[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)`
